# ---------------------------------------------------------------------------
# 01b - Gap detection and repair
#
# A completeness audit on the downloaded data found seven local dates on which
# only one of the two VIIRS platforms contributed any row across the whole of
# Kalimantan. A platform observing an island-sized box and returning nothing is
# not a physical result; it is a request that was served empty. On those dates the
# detection count is halved by construction, which would appear in the daily
# series as a decline in fire activity.
#
# This script finds such dates, re-queries only the affected platform and dates,
# merges the result, and writes a report of what changed. It is separated from
# R/01 so the repair is visible in the pipeline rather than buried in a re-run.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(glue)
})
source("R/00_config.R")

map_key <- Sys.getenv("FIRMS_MAP_KEY", "")
stopifnot(nzchar(map_key))
KBB <- unname(KALIMANTAN_BBOX[c("west","south","east","north")])

firms_get <- function(source, days, start, tries = 4L) {
  url <- glue("https://firms.modaps.eosdis.nasa.gov/api/area/csv/{map_key}/{source}/",
              "{paste(KBB, collapse=',')}/{days}/{format(start, '%Y-%m-%d')}")
  for (k in seq_len(tries)) {
    out <- try(suppressWarnings(read_csv(url, show_col_types = FALSE, progress = FALSE)),
               silent = TRUE)
    if (!inherits(out, "try-error") && "latitude" %in% names(out) && nrow(out) > 0)
      return(out)
    Sys.sleep(2 * k)
  }
  tibble()
}

raw <- readRDS(file.path(DIR_DERIVED, "firms_raw.rds"))
n_before <- nrow(raw)

# --- detect: which (date, platform) pairs are absent? ----------------------
all_days <- seq(WINDOW_START, WINDOW_END, by = "day")
present <- raw |> mutate(dt = as.Date(acq_date)) |> distinct(dt, platform)
expected <- expand.grid(dt = all_days, platform = c("S-NPP", "NOAA-20"),
                        stringsAsFactors = FALSE) |> as_tibble()
missing <- anti_join(expected, present, by = c("dt", "platform")) |> arrange(dt, platform)
message(sprintf("(date, platform) combinations absent: %d", nrow(missing)))
if (nrow(missing)) print(as.data.frame(missing))

if (!nrow(missing)) {
  message("no gaps; nothing to repair")
} else {
  src <- c("S-NPP" = "VIIRS_SNPP_NRT", "NOAA-20" = "VIIRS_NOAA20_NRT")
  # Group consecutive missing dates per platform into blocks of at most 5 days,
  # the API's documented maximum, so the repair costs as few requests as possible.
  blocks <- missing |> group_by(platform) |>
    mutate(grp = cumsum(c(1, diff(dt) != 1))) |>
    group_by(platform, grp) |>
    summarise(start = min(dt), days = as.integer(max(dt) - min(dt)) + 1L,
              .groups = "drop") |>
    mutate(days = pmin(days, FIRMS_MAX_DAY_RANGE))
  print(as.data.frame(blocks))

  filled <- pmap_dfr(list(blocks$platform, blocks$start, blocks$days),
                     function(pf, st, dy) {
    d <- firms_get(src[[pf]], dy, st)
    message(glue("  refetch {pf} {format(st, '%d %b')} +{dy}d: {nrow(d)} rows"))
    if (!nrow(d)) return(NULL)
    d |> mutate(platform = pf, acq_source = "NRT")
  })

  if (nrow(filled)) {
    raw <- bind_rows(raw, filled) |>
      distinct(platform, acq_date, acq_time, latitude, longitude, .keep_all = TRUE)
  }
  message(sprintf("rows: %s -> %s (added %s)", format(n_before, big.mark = ","),
                  format(nrow(raw), big.mark = ","),
                  format(nrow(raw) - n_before, big.mark = ",")))
  saveRDS(raw, file.path(DIR_DERIVED, "firms_raw.rds"))
}

# --- verify: re-audit after repair -----------------------------------------
present2 <- raw |> mutate(dt = as.Date(acq_date)) |> distinct(dt, platform)
still <- anti_join(expected, present2, by = c("dt", "platform")) |> arrange(dt, platform)
per_day <- present2 |> count(dt, name = "n_platforms")
report <- list(
  n_before = n_before, n_after = nrow(raw),
  missing_before = missing, missing_after = still,
  days_with_one_platform_before = nrow(missing),
  days_with_one_platform_after = nrow(still),
  per_day = per_day)
saveRDS(report, file.path(DIR_DERIVED, "gap_report.rds"))
message(sprintf("after repair, (date, platform) combinations still absent: %d", nrow(still)))
if (nrow(still)) print(as.data.frame(still))
