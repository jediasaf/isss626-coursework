# ---------------------------------------------------------------------------
# 01_prepare.R - import, audit, clean, join and aggregate
#
# Output: data/derived/prepared.rds, a list holding
#   tsp        township polygons (sf, projected) with total counts per variable
#   ev         the retained political-violence events (sf, projected)
#   audit      a step-by-step record of every row removed and why
#   names_chk  ACLED admin3 names versus the township each point falls in
#   cube       township x quarter counts (complete, zeros filled)
#   cube_m     township x month counts, for the sensitivity check
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr); library(lubridate)
})
source("R/00_config.R")
sf_use_s2(FALSE)

## 1. Boundaries --------------------------------------------------------------
if (!file.exists(MIMU_FILE)) {
  dir.create(dirname(MIMU_FILE), recursive = TRUE, showWarnings = FALSE)
  zip <- file.path(dirname(MIMU_FILE), "mimu_adm3.zip")
  download.file(MIMU_WFS, zip, mode = "wb")
  unzip(zip, exdir = dirname(MIMU_FILE))
}
tsp_raw <- st_read(MIMU_FILE, quiet = TRUE)
stopifnot(all(tsp_raw$PCode_V == "9.4"), !anyDuplicated(tsp_raw$TS_PCODE))

# Sorting by PCode matters later: sfdep 0.2.5 emerging_hotspot_analysis() pastes
# its results against the geometry order after splitting by location
# alphabetically. With the geometry already in PCode order the two agree.
tsp <- tsp_raw |>
  select(ST, ST_PCODE, DT, DT_PCODE, TS, TS_PCODE) |>
  arrange(TS_PCODE) |>
  st_make_valid() |>
  st_transform(STUDY_CRS)
tsp$area_km2 <- as.numeric(st_area(tsp)) / 1e6

## 2. Events -------------------------------------------------------------------
raw <- read_csv(ACLED_FILE, show_col_types = FALSE,
                col_types = cols(.default = col_character()))
raw <- raw |>
  mutate(event_date = as.Date(event_date),
         across(c(latitude, longitude), as.numeric),
         across(c(geo_precision, time_precision, fatalities, year), as.integer))

audit <- list()
step <- function(df, label, why) {
  audit[[length(audit) + 1]] <<- tibble(step = label, rows = nrow(df), why = why)
  df
}

ev <- raw |> step("ACLED extract as supplied", "every row in the eLearn file")

stopifnot(!anyDuplicated(ev$event_id_cnty))
ev <- ev |>
  filter(event_date >= WINDOW_START, event_date <= WINDOW_END) |>
  step("Inside the study window", "1 Jan 2021 to 30 Sep 2025, checked rather than assumed")

ev <- ev |>
  filter(event_type %in% PV_TYPES) |>
  step("Political violence only",
       "battles, explosions and remote violence, violence against civilians")

ev <- ev |>
  filter(country == "Myanmar") |>
  step("Country is Myanmar", "removes events ACLED codes to the Indian Ocean")

# Kept for the audit: how many events lacked a usable location.
ev <- ev |>
  filter(geo_precision %in% GEO_PRECISION_KEEP) |>
  step("Location known to township or better",
       "geo_precision 3 events sit on a regional centroid")

## 3. Point-in-polygon join ---------------------------------------------------
ev_sf <- st_as_sf(ev, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
  st_transform(STUDY_CRS)
idx <- st_intersects(ev_sf, tsp)
n_hit <- lengths(idx)
# A point on a shared boundary intersects two townships; take the first so it
# is counted once.
ev_sf$tsp_i <- vapply(idx, function(v) if (length(v)) v[1] else NA_integer_, 1L)
on_boundary <- sum(n_hit > 1)

miss <- which(is.na(ev_sf$tsp_i))
near <- st_nearest_feature(ev_sf[miss, ], tsp)
d_near <- as.numeric(st_distance(ev_sf[miss, ], tsp[near, ], by_element = TRUE))
snap <- d_near <= SNAP_TOLERANCE_M
ev_sf$tsp_i[miss[snap]] <- near[snap]
snapped <- tibble(event_id_cnty = ev_sf$event_id_cnty[miss], admin3 = ev_sf$admin3[miss],
                  location = ev_sf$location[miss], dist_km = d_near / 1000, snapped = snap)

ev_sf <- ev_sf |> filter(!is.na(tsp_i)) |>
  step("Inside a township, or within 5 km of one",
       "points beyond the 1:250k coastline are offshore")
ev_sf$TS_PCODE <- tsp$TS_PCODE[ev_sf$tsp_i]
ev_sf$TS       <- tsp$TS[ev_sf$tsp_i]
ev_sf$ST       <- tsp$ST[ev_sf$tsp_i]

## 4. Name agreement as a diagnostic, not a join key ---------------------------
# ACLED records its own admin3 name. Joining on names would fail wherever the
# two sources romanise differently, and would silently lose events. The spatial
# join is used instead, and the names are compared afterwards to measure how
# much a name join would have got wrong.
norm <- function(x) tolower(gsub("[^A-Za-z]", "", x))
names_chk <- st_drop_geometry(ev_sf) |>
  mutate(match = norm(admin3) == norm(TS)) |>
  count(admin3, TS, match, name = "events") |>
  arrange(desc(events))

## 5. Variables ----------------------------------------------------------------
ev_sf <- ev_sf |>
  mutate(var = case_when(event_type == "Battles" ~ "battles",
                         event_type == "Explosions/Remote violence" ~ "erv",
                         event_type == "Violence against civilians" ~ "vac"),
         quarter = floor_date(event_date, "quarter"),
         month   = floor_date(event_date, "month"))

totals <- st_drop_geometry(ev_sf) |>
  group_by(TS_PCODE) |>
  summarise(pv = n(), battles = sum(var == "battles"), erv = sum(var == "erv"),
            vac = sum(var == "vac"), fatalities = sum(fatalities), .groups = "drop")
tsp <- tsp |> left_join(totals, by = "TS_PCODE") |>
  mutate(across(c(pv, battles, erv, vac, fatalities), \(x) coalesce(x, 0L)))

## 6. Space-time cubes ----------------------------------------------------------
# complete() fills every township-period that recorded nothing with an explicit
# zero. A cube with missing cells is not a cube, and sfdep refuses it.
make_cube <- function(period_col, periods) {
  st_drop_geometry(ev_sf) |>
    rename(period = all_of(period_col)) |>
    group_by(TS_PCODE, period) |>
    summarise(pv = n(), battles = sum(var == "battles"), erv = sum(var == "erv"),
              vac = sum(var == "vac"), fatalities = sum(fatalities), .groups = "drop") |>
    complete(TS_PCODE = tsp$TS_PCODE, period = periods,
             fill = list(pv = 0L, battles = 0L, erv = 0L, vac = 0L, fatalities = 0L)) |>
    arrange(period, TS_PCODE)
}
quarters <- seq(floor_date(WINDOW_START, "quarter"), floor_date(WINDOW_END, "quarter"), by = "quarter")
months   <- seq(WINDOW_START, floor_date(WINDOW_END, "month"), by = "month")
cube   <- make_cube("quarter", quarters)
cube_m <- make_cube("month", months)
stopifnot(nrow(cube) == nrow(tsp) * length(quarters),
          sum(cube$pv) == nrow(ev_sf), sum(cube_m$pv) == nrow(ev_sf))

## 7. Diagnostics carried into the report -------------------------------------
qa <- list(
  n_raw = nrow(raw),
  by_type_raw = count(raw, event_type),
  geo3_pv = sum(raw$event_type %in% PV_TYPES & raw$geo_precision == 3L & raw$country == "Myanmar"),
  indian_ocean = sum(raw$country != "Myanmar"),
  admin3_na_pv = sum(is.na(raw$admin3) & raw$event_type %in% PV_TYPES),
  time_prec = count(filter(raw, event_type %in% PV_TYPES), time_precision),
  on_boundary = on_boundary,
  snapped = snapped,
  fat_top = st_drop_geometry(ev_sf) |> arrange(desc(fatalities)) |>
    select(event_id_cnty, event_date, sub_event_type, TS, fatalities) |> head(10),
  fat_share_top1pct = {
    f <- sort(ev_sf$fatalities, decreasing = TRUE)
    sum(f[seq_len(ceiling(length(f) * 0.01))]) / sum(f)
  },
  zero_tsp = sum(tsp$pv == 0),
  monthly = st_drop_geometry(ev_sf) |> count(month, var),
  sub_types = st_drop_geometry(ev_sf) |> count(var, sub_event_type, sort = TRUE),
  name_mismatch_events = sum(names_chk$events[!names_chk$match & !is.na(names_chk$match)]),
  name_na_events = sum(names_chk$events[is.na(names_chk$match)])
)

audit <- bind_rows(audit) |> mutate(removed = lag(rows) - rows)

ev_keep <- ev_sf |> select(event_id_cnty, event_date, quarter, month, event_type,
                           sub_event_type, var, actor1, actor2, fatalities,
                           geo_precision, time_precision, TS_PCODE, TS, ST)

saveRDS(list(tsp = tsp, ev = ev_keep, audit = audit, names_chk = names_chk,
             cube = cube, cube_m = cube_m, quarters = quarters, months = months, qa = qa),
        file.path(DIR_DERIVED, "prepared.rds"))

print(audit)
cat("on boundary:", on_boundary, " snapped:", sum(snap), " dropped offshore:", sum(!snap), "\n")
cat("name mismatch events:", qa$name_mismatch_events, " admin3 NA events:", qa$name_na_events, "\n")
cat("townships with zero events:", qa$zero_tsp, "\n")
