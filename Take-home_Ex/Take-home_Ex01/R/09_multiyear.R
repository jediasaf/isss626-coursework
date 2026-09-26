# ---------------------------------------------------------------------------
# 09 - Multi-year context
#
# The single largest limitation of the core analysis was that it covered 87 days
# of one season, so it could not say whether 2026 was unusual. I had recorded
# that as a hard constraint, on the evidence that the FIRMS standard-processing
# (SP) product returned zero rows for 2026. That inference was wrong: SP is
# empty for 2026 only because near-real-time granules are replaced by standard
# science quality after a lag of about five months. For earlier years the SAME
# `area` REST endpoint serves the SP archive normally, back to the start of the
# VIIRS record.
#
# So the multi-year record was scriptable from the same key the whole time. This
# script pulls it. Two products:
#
#   A  the identical 1 July to 25 September window for every year from 2013,
#      which is what makes 2026 comparable rather than merely described
#   B  two complete calendar years, to establish the wet-season baseline the
#      87-day window cannot contain
#
# S-NPP only. It is the one platform continuous across the whole period, so
# counts are comparable year to year; adding NOAA-20 from 2018 would introduce a
# step change in detection capability that would be indistinguishable from a
# change in fire activity.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(glue); library(sf)
})
source("R/00_config.R")
set.seed(SEED)

map_key <- Sys.getenv("FIRMS_MAP_KEY", "")
stopifnot(nzchar(map_key))
KBB <- unname(KALIMANTAN_BBOX[c("west","south","east","north")])

# Pulang Pisau bounding box, to keep the pull small. The clip to the polygon
# happens afterwards, exactly as in R/02.
win_sf <- readRDS(file.path(DIR_DERIVED, "kal_adm2.rds")) |>
  filter(shapeName == STUDY_REGENCY)
bb <- st_bbox(win_sf)
PBB <- c(bb[["xmin"]] - 0.05, bb[["ymin"]] - 0.05, bb[["xmax"]] + 0.05, bb[["ymax"]] + 0.05)

firms_get <- function(source, bbox, days, start, tries = 3L) {
  url <- glue("https://firms.modaps.eosdis.nasa.gov/api/area/csv/{map_key}/{source}/",
              "{paste(bbox, collapse=',')}/{days}/{format(start, '%Y-%m-%d')}")
  for (k in seq_len(tries)) {
    # Read every column as character. Across 300-odd requests a few chunks come
    # back with a column readr guesses differently (latitude as character, for
    # instance), and bind_rows then refuses to combine them. Coercing once at
    # the end is safer than trusting per-chunk type inference.
    out <- try(suppressWarnings(read_csv(url, show_col_types = FALSE, progress = FALSE,
                                         col_types = cols(.default = col_character()))),
               silent = TRUE)
    if (!inherits(out, "try-error") && "latitude" %in% names(out)) {
      return(structure(out, served = TRUE))
    }
    Sys.sleep(2 * k)
  }
  structure(tibble(), served = FALSE)
}

pull_range <- function(from, to, source, label) {
  starts <- seq(from, to, by = paste(FIRMS_MAX_DAY_RANGE, "days"))
  map_dfr(starts, function(s) {
    days <- as.integer(min(to, s + FIRMS_MAX_DAY_RANGE - 1L) - s) + 1L
    d <- firms_get(source, PBB, days, s)
    if (!isTRUE(attr(d, "served"))) {
      message(glue("  !! {label} {s}: NOT SERVED"))
      return(tibble(block_start = s, served = FALSE, n = NA_integer_))
    }
    if (!nrow(d)) return(NULL)
    # Keep only what the analysis needs, so 300 chunks stay small in memory.
    keep <- intersect(c("latitude","longitude","acq_date","acq_time","confidence",
                        "frp","daynight","satellite"), names(d))
    d <- d[, keep, drop = FALSE]
    d$block_start <- as.character(s)
    d
  })
}

# One coercion, applied after everything is bound.
coerce_firms <- function(d) {
  if (!nrow(d)) return(d)
  d |> mutate(
    latitude  = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude)),
    frp       = suppressWarnings(as.numeric(frp)),
    acq_date  = as.Date(acq_date))
}

# --- A. the same 87-day window, every year ---------------------------------
YEARS <- 2013:2025
message("A. pulling 1 Jul to 25 Sep for each year 2013-2025 (S-NPP, SP archive)")
annual <- map_dfr(YEARS, function(y) {
  from <- as.Date(sprintf("%d-07-01", y)); to <- as.Date(sprintf("%d-09-25", y))
  d <- pull_range(from, to, "VIIRS_SNPP_SP", glue("{y}"))
  if (!nrow(d) || !"latitude" %in% names(d)) {
    message(glue("  {y}: 0 rows")); return(tibble(year = y, n_raw = 0L))
  }
  d$year <- y
  message(glue("  {y}: {nrow(d)} rows"))
  d
})
annual <- coerce_firms(annual)
saveRDS(annual, file.path(DIR_DERIVED, "multiyear_window_raw.rds"))

# --- B. two complete years, for the seasonal baseline -----------------------
message("B. pulling complete calendar years 2023 and 2024 (S-NPP, SP archive)")
fullyear <- map_dfr(c(2023, 2024), function(y) {
  d <- pull_range(as.Date(sprintf("%d-01-01", y)), as.Date(sprintf("%d-12-31", y)),
                  "VIIRS_SNPP_SP", glue("full{y}"))
  if (!nrow(d) || !"latitude" %in% names(d)) return(tibble(year = y, n_raw = 0L))
  d$year <- y
  message(glue("  {y}: {nrow(d)} rows"))
  d
})
fullyear <- coerce_firms(fullyear)
saveRDS(fullyear, file.path(DIR_DERIVED, "multiyear_fullyear_raw.rds"))
message("multi-year acquisition complete")
