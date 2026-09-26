# ---------------------------------------------------------------------------
# 10 - Is 2026 unusual, and what does a wet season look like?
#
# Two comparability rules, both of which cost something and both of which are
# worth the cost:
#
#  1. S-NPP only, every year including 2026. NOAA-20 joins the record in 2018,
#     so a two-platform count for recent years against a one-platform count for
#     earlier years would show a step change in detection capability that is
#     indistinguishable from a change in fire activity. This means the 2026
#     figure here is SMALLER than the headline figure in the core analysis,
#     which uses both platforms.
#
#  2. The identical calendar window, 1 July to 25 September, for every year.
#     Comparing a full season against an 87-day slice would confound length
#     with intensity.
#
# The caveat I cannot remove: 2026 comes from the near-real-time product and
# every earlier year comes from the reprocessed standard-quality product. The
# two differ slightly. I cannot quantify the difference, because the retention
# windows of the two products do not overlap anywhere.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(purrr);
  library(dplyr); library(sf); library(glue); library(tidyr)
})
source("R/00_config.R")
set.seed(SEED)

P <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
win_p <- P$win_p

# Same screening as R/02: confidence, coordinates, clip, then collapse to the
# 375 m cell so the unit matches the core analysis.
prep_year <- function(d) {
  if (!nrow(d)) return(tibble())
  d |>
    filter(!is.na(latitude), !is.na(longitude),
           between(latitude, -90, 90), between(longitude, -180, 180),
           tolower(substr(confidence, 1, 1)) %in% c("n", "h")) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    st_transform(STUDY_CRS) |>
    (\(x) x[lengths(st_within(x, win_p)) > 0, ])() |>
    (\(x) {
      xy <- st_coordinates(x)
      x |> st_drop_geometry() |>
        mutate(cell = paste(floor(xy[,1] / DEDUP_GRID_M),
                            floor(xy[,2] / DEDUP_GRID_M), sep = "_"))
    })()
}

# --- A. the same window, every year ----------------------------------------
annual_raw <- readRDS(file.path(DIR_DERIVED, "multiyear_window_raw.rds"))
annual <- annual_raw |> group_split(year) |>
  map_dfr(function(d) {
    y <- d$year[1]
    p <- prep_year(d)
    if (!nrow(p)) return(tibble(year = y, detections = 0L, cells = 0L, cell_days = 0L))
    tibble(year = y,
           detections = nrow(p),
           cells      = n_distinct(p$cell),
           cell_days  = nrow(distinct(p, cell, acq_date)))
  })

# 2026, recomputed S-NPP only so it is comparable with the history above.
f26 <- P$fires |> st_drop_geometry() |> filter(platform == "S-NPP")
annual <- bind_rows(annual, tibble(
  year = 2026L,
  detections = nrow(f26),
  cells = n_distinct(f26$cell),
  cell_days = nrow(distinct(f26, cell, date_local))
)) |> arrange(year)

annual <- annual |>
  mutate(rank_cells = rank(-cells, ties.method = "min"),
         vs_median  = cells / median(cells[year != 2026]))
print(as.data.frame(annual))

# --- B. the seasonal cycle, from complete years ----------------------------
fy_raw <- readRDS(file.path(DIR_DERIVED, "multiyear_fullyear_raw.rds"))
seasonal <- fy_raw |> group_split(year) |>
  map_dfr(function(d) {
    y <- d$year[1]; p <- prep_year(d)
    if (!nrow(p)) return(tibble())
    p |> mutate(month = as.integer(format(acq_date, "%m"))) |>
      distinct(cell, acq_date, month) |>
      count(year = y, month, name = "cell_days")
  }) |>
  complete(year, month = 1:12, fill = list(cell_days = 0L))
print(as.data.frame(seasonal |> pivot_wider(names_from = month, values_from = cell_days)))

# The quantity the 87-day window cannot contain: how much of a year's fire
# happens inside it, and what the wet-season floor actually is.
amplitude <- seasonal |> group_by(year) |>
  summarise(annual_total = sum(cell_days),
            in_window_JulSep = sum(cell_days[month %in% 7:9]),
            share_JulSep = in_window_JulSep / annual_total,
            wet_JanApr = sum(cell_days[month %in% 1:4]),
            peak_month = month[which.max(cell_days)],
            peak_over_wet = max(cell_days) / pmax(mean(cell_days[month %in% 1:4]), 0.5),
            .groups = "drop")
print(as.data.frame(amplitude))

saveRDS(list(annual = annual, seasonal = seasonal, amplitude = amplitude,
             years_pulled = sort(unique(annual_raw$year))),
        file.path(DIR_DERIVED, "multiyear.rds"))
message("multi-year analysis cached")
