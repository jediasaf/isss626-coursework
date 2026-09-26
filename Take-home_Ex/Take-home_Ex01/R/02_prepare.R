# ---------------------------------------------------------------------------
# 02 - Preparation, quality assessment, and construction of the point patterns
#
# Every exclusion is counted and written to data/derived/qa_steps.rds so the
# report can show an audit trail rather than assert that cleaning happened.
#
# Three patterns are built, because "one VIIRS detection" is not "one fire":
#   ppp_all      every retained detection            (sensitivity only)
#   ppp_daycell  one per 375 m cell per day          (spatio-temporal unit)
#   ppp_cell     one per 375 m cell over the window  (primary spatial unit)
# The rationale is in the report; the short version is that a single large,
# multi-day fire yields many detections, so the raw pattern measures sensor
# sampling as much as fire occurrence.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr); library(glue)
  library(spatstat.geom); library(units)
})
source("R/00_config.R")
set.seed(SEED)

qa <- tibble(step = character(), rule = character(),
             n_in = integer(), n_removed = integer(), n_out = integer())
record <- function(step, rule, n_in, n_out) {
  qa <<- add_row(qa, step = step, rule = rule, n_in = as.integer(n_in),
                 n_removed = as.integer(n_in - n_out), n_out = as.integer(n_out))
  message(glue("{step}: {n_in} -> {n_out} (removed {n_in - n_out})"))
}

raw <- readRDS(file.path(DIR_DERIVED, "firms_raw.rds"))
kal <- readRDS(file.path(DIR_DERIVED, "kal_adm2.rds"))
n0 <- nrow(raw)
record("00 downloaded", "rows returned by the FIRMS area API", n0, n0)

# --- 1. field typing and normalisation -------------------------------------
# confidence: the area API returns l/n/h, the regional files return words.
norm_conf <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    substr(x, 1, 1) == "l" ~ "low",
    substr(x, 1, 1) == "n" ~ "nominal",
    substr(x, 1, 1) == "h" ~ "high",
    TRUE ~ NA_character_)
}
d <- raw |>
  mutate(
    acq_date   = as.Date(acq_date),
    # acq_time is HHMM without a leading zero; pad before parsing.
    acq_hhmm   = sprintf("%04d", suppressWarnings(as.integer(acq_time))),
    acq_hour   = as.integer(substr(acq_hhmm, 1, 2)),
    acq_min    = as.integer(substr(acq_hhmm, 3, 4)),
    # FIRMS timestamps are UTC. Kalimantan runs WITA, UTC+8. Local time is what
    # matters for interpreting the diurnal cycle and for the definition of a
    # "fire day", so both are carried.
    ts_utc     = as.POSIXct(paste(acq_date, acq_hhmm), format = "%Y-%m-%d %H%M", tz = "UTC"),
    ts_local   = ts_utc + 8 * 3600,
    date_local = as.Date(ts_local),
    confidence = norm_conf(confidence),
    daynight   = ifelse(toupper(substr(daynight, 1, 1)) == "D", "day", "night"),
    frp        = suppressWarnings(as.numeric(frp))
  )
has_type <- "type" %in% names(d)
message(glue("`type` column present in this product: {has_type}"))

# --- 2. coordinate validity -------------------------------------------------
n_in <- nrow(d)
d <- d |> filter(
  !is.na(latitude), !is.na(longitude),
  is.finite(latitude), is.finite(longitude),
  between(latitude, -90, 90), between(longitude, -180, 180),
  !(latitude == 0 & longitude == 0)          # null-island sentinel
)
record("01 coordinates", "non-missing, finite, in range, not (0,0)", n_in, nrow(d))

# --- 3. timestamp validity --------------------------------------------------
n_in <- nrow(d)
d <- d |> filter(!is.na(ts_utc), !is.na(acq_hour), between(acq_hour, 0, 23),
                 between(acq_min, 0, 59))
record("02 timestamps", "parseable UTC timestamp, valid HH and MM", n_in, nrow(d))

# --- 4. confidence screen ---------------------------------------------------
conf_tab <- d |> count(confidence)
n_in <- nrow(d)
d <- d |> filter(confidence %in% CONF_KEEP)
record("03 confidence", "drop low-confidence detections", n_in, nrow(d))

# --- 5. source type, where the product carries it ---------------------------
if (has_type) {
  n_in <- nrow(d)
  d <- d |> filter(type %in% TYPE_KEEP)
  record("04 source type", "type == 0 (presumed vegetation fire)", n_in, nrow(d))
} else {
  record("04 source type", "column absent from NRT product; see persistence screen",
         nrow(d), nrow(d))
}

# --- 6. temporal window -----------------------------------------------------
n_in <- nrow(d)
d <- d |> filter(date_local >= WINDOW_START, date_local <= WINDOW_END)
record("05 temporal window", glue("{WINDOW_START} to {WINDOW_END}, local date"), n_in, nrow(d))

# --- 7. exact duplicate records --------------------------------------------
# The same granule can be delivered twice. A duplicate is the same platform,
# the same instant, and the same coordinates to full reported precision.
n_in <- nrow(d)
d <- d |> distinct(platform, ts_utc, latitude, longitude, .keep_all = TRUE)
record("06 exact duplicates", "unique (platform, timestamp, lat, lon)", n_in, nrow(d))

# --- 8. spatial clip to the study window -----------------------------------
win_sf <- kal |> filter(shapeName == STUDY_REGENCY) |> st_make_valid()
stopifnot(nrow(win_sf) == 1, all(st_is_valid(win_sf)))
win_p  <- st_transform(win_sf, STUDY_CRS)

pts <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
  st_transform(STUDY_CRS)
n_in <- nrow(pts)
inside <- lengths(st_within(pts, win_p)) > 0
pts <- pts[inside, ]
record("07 spatial clip", glue("within the {STUDY_REGENCY} polygon (also removes offshore)"),
       n_in, nrow(pts))

# --- 9. persistence screen, standing in for the absent `type` field --------
# Assign every detection to a 375 m cell of a fixed UTM grid, then count the
# number of distinct local dates on which each cell was detected. A cell lit on
# most days of an 87-day window is a fixed installation, not a vegetation fire.
xy <- st_coordinates(pts)
pts <- pts |> mutate(
  cell_x = floor(xy[, 1] / DEDUP_GRID_M),
  cell_y = floor(xy[, 2] / DEDUP_GRID_M),
  cell   = paste(cell_x, cell_y, sep = "_")
)
n_days_window <- as.integer(WINDOW_END - WINDOW_START) + 1L
cell_days <- pts |> st_drop_geometry() |>
  group_by(cell) |>
  summarise(n_det = n(), n_days = n_distinct(date_local),
            mean_frp = mean(frp, na.rm = TRUE), .groups = "drop") |>
  mutate(day_frac = n_days / n_days_window)

# Sensitivity across candidate thresholds. The grid deliberately extends well
# below the adopted threshold, so the table shows where the screen would begin
# to bite rather than only that it is quiet at the chosen value. A screen that
# removes nothing is only informative alongside evidence that it is functional.
persist_sens <- tibble(threshold = c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 0.80)) |>
  rowwise() |>
  mutate(n_cells = sum(cell_days$day_frac > threshold),
         n_detections = sum(cell_days$n_det[cell_days$day_frac > threshold])) |>
  ungroup()

# Second, independent check on the same question. A fixed installation such as a
# flare stack emits a high and steady radiative power; smouldering peat emits a
# low one. If the most persistent cells were industrial, their mean FRP would be
# an order of magnitude above the window median.
persist_top <- cell_days |> arrange(desc(n_days)) |> head(10) |>
  transmute(cell, n_det, n_days, day_frac = round(day_frac, 3),
            mean_frp = round(mean_frp, 1))
frp_median_all <- median(pts$frp, na.rm = TRUE)
persist_max_days <- max(cell_days$n_days)

persistent_cells <- cell_days$cell[cell_days$day_frac > PERSIST_DAY_FRAC]
n_in <- nrow(pts)
pts_flagged <- pts |> filter(cell %in% persistent_cells)
pts <- pts |> filter(!cell %in% persistent_cells)
record("08 persistence screen",
       glue("drop cells detected on >{round(PERSIST_DAY_FRAC*100)}% of {n_days_window} days"),
       n_in, nrow(pts))

# --- 10. the three analysis patterns ---------------------------------------
win <- as.owin(st_geometry(win_p))
Window_area_km2 <- as.numeric(set_units(st_area(win_p), "km^2"))

# Representative detection per cell-day and per cell: the earliest observation,
# a deterministic rule that preserves a genuine sensor coordinate rather than
# snapping to a lattice, which would manufacture regularity at small distances.
daycell <- pts |> arrange(ts_utc) |> distinct(cell, date_local, .keep_all = TRUE)
cellone <- pts |> arrange(ts_utc) |> distinct(cell, .keep_all = TRUE)

mk_ppp <- function(x, label) {
  co <- st_coordinates(x)
  p <- ppp(co[, 1], co[, 2], window = win, checkdup = FALSE)
  nrej <- sum(!inside.owin(co[, 1], co[, 2], win))
  message(glue("{label}: n = {npoints(p)}, rejected on the boundary = {nrej}, ",
               "duplicated coordinate pairs = {sum(duplicated(co))}"))
  marks(p) <- x |> st_drop_geometry() |>
    transmute(frp, confidence = factor(confidence, levels = c("nominal","high")),
              daynight = factor(daynight), platform = factor(platform),
              date_local, doy = as.integer(date_local - WINDOW_START) + 1L)
  p
}
ppp_all     <- mk_ppp(pts,     "ppp_all")
ppp_daycell <- mk_ppp(daycell, "ppp_daycell")
ppp_cell    <- mk_ppp(cellone, "ppp_cell")

# --- 11. diagnostics used in the report ------------------------------------
multiplicity_tab <- cell_days |>
  filter(!cell %in% persistent_cells) |>
  count(n_det, name = "n_cells") |> arrange(n_det)

daily <- pts |> st_drop_geometry() |>
  count(date_local, name = "n_detections") |>
  right_join(tibble(date_local = seq(WINDOW_START, WINDOW_END, by = "day")),
             by = "date_local") |>
  mutate(n_detections = tidyr::replace_na(n_detections, 0L)) |>
  arrange(date_local)
daily_daycell <- daycell |> st_drop_geometry() |>
  count(date_local, name = "n_events") |>
  right_join(tibble(date_local = seq(WINDOW_START, WINDOW_END, by = "day")),
             by = "date_local") |>
  mutate(n_events = tidyr::replace_na(n_events, 0L)) |> arrange(date_local)

out <- list(
  fires = pts, flagged = pts_flagged, win_sf = win_sf, win_p = win_p, win = win,
  ppp_all = ppp_all, ppp_daycell = ppp_daycell, ppp_cell = ppp_cell,
  qa = qa, conf_tab = conf_tab, cell_days = cell_days,
  persist_sens = persist_sens, persistent_cells = persistent_cells,
  persist_top = persist_top, frp_median_all = frp_median_all,
  persist_max_days = persist_max_days,
  multiplicity = multiplicity_tab, daily = daily, daily_daycell = daily_daycell,
  area_km2 = Window_area_km2, has_type = has_type, n_days = n_days_window
)
saveRDS(out, file.path(DIR_DERIVED, "prepared.rds"))
saveRDS(qa,  file.path(DIR_DERIVED, "qa_steps.rds"))

message("\n--- audit trail ---"); print(as.data.frame(qa), row.names = FALSE)
message(glue("\nwindow area: {round(Window_area_km2)} km2"))
message(glue("ppp_all {npoints(ppp_all)} | ppp_daycell {npoints(ppp_daycell)} | ppp_cell {npoints(ppp_cell)}"))
message(glue("mean intensity, ppp_cell: {signif(npoints(ppp_cell)/Window_area_km2, 3)} per km2"))
