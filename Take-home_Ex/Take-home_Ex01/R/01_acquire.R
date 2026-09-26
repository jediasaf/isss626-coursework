# ---------------------------------------------------------------------------
# 01 - Data acquisition
#
# Sources, both re-creatable by running this file:
#   (a) geoBoundaries gbOpen Indonesia ADM1 / ADM2, CC BY 3.0 IGO
#   (b) NASA FIRMS VIIRS 375 m Collection 2 active-fire detections, area API
#
# FIRMS needs a free MAP_KEY (https://firms.modaps.eosdis.nasa.gov/api/map_key/).
# Put it in .Renviron as FIRMS_MAP_KEY. The key is deliberately not committed.
#
# Outputs (data/derived/):
#   kal_adm2.rds, kal_prov.rds, firms_raw.rds,
#   regency_ranking.rds, retention_probe.rds, acquisition_log.md
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(purrr); library(glue)
  library(units)
})
source("R/00_config.R")
dir.create(DIR_RAW, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_DERIVED, showWarnings = FALSE, recursive = TRUE)
sf_use_s2(TRUE)

log_lines <- character(0)
say <- function(..., .envir = parent.frame()) {
  msg <- as.character(glue(..., .envir = .envir))
  log_lines <<- c(log_lines, msg); message(msg); invisible(msg)
}

map_key <- Sys.getenv("FIRMS_MAP_KEY", "")
if (!nzchar(map_key)) stop(
  "FIRMS_MAP_KEY is not set. Copy .Renviron.example to .Renviron and add a free ",
  "key from https://firms.modaps.eosdis.nasa.gov/api/map_key/")

# --- (a) Administrative boundaries -----------------------------------------
GB <- list(
  ADM1 = "https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen/IDN/ADM1/geoBoundaries-IDN-ADM1.geojson",
  ADM2 = "https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen/IDN/ADM2/geoBoundaries-IDN-ADM2.geojson"
)
for (lvl in names(GB)) {
  f <- file.path(DIR_RAW, glue("geoBoundaries-IDN-{lvl}.geojson"))
  if (!file.exists(f)) {
    say("downloading geoBoundaries IDN {lvl}")
    download.file(GB[[lvl]], f, mode = "wb", quiet = TRUE)
  }
  say("geoBoundaries IDN {lvl}: {round(file.size(f)/1e6,1)} MB, pinned commit 9469f09")
}

adm1 <- st_read(file.path(DIR_RAW, "geoBoundaries-IDN-ADM1.geojson"), quiet = TRUE)
kal_prov <- adm1 |> filter(grepl("Kalimantan", shapeName, ignore.case = TRUE)) |> st_make_valid()
say("Kalimantan provinces: {paste(sort(kal_prov$shapeName), collapse=', ')}")

# Read only the regencies intersecting the Kalimantan bounding box, so the
# 159 MB national file is never held in memory in full.
adm2 <- st_read(file.path(DIR_RAW, "geoBoundaries-IDN-ADM2.geojson"),
                wkt_filter = st_as_text(st_as_sfc(st_bbox(kal_prov))), quiet = TRUE) |>
  st_make_valid()

# Retain a regency if a representative interior point falls inside the
# provincial union. Projecting first removes the geographic-coordinate warning
# from st_point_on_surface and makes the test planar and unambiguous.
kal_union_p <- st_union(st_transform(kal_prov, STUDY_CRS))
adm2_p      <- st_transform(adm2, STUDY_CRS)
reps        <- st_point_on_surface(st_geometry(adm2_p))
kal_adm2    <- adm2[lengths(st_intersects(reps, kal_union_p)) > 0, ]
say("Kalimantan regencies/cities retained: {nrow(kal_adm2)} of {nrow(adm2)} in bbox")

# Parent province by containment of the representative point (not nearest
# feature, which would mis-assign a regency lying just outside a province edge).
kal_reps <- st_point_on_surface(st_geometry(st_transform(kal_adm2, STUDY_CRS)))
kal_prov_p <- st_transform(kal_prov, STUDY_CRS)
hit <- st_within(kal_reps, kal_prov_p, sparse = TRUE)
kal_adm2$province <- vapply(hit, function(i)
  if (length(i)) kal_prov$shapeName[i[1]] else NA_character_, character(1))
say("regencies with an unresolved parent province: {sum(is.na(kal_adm2$province))}")

saveRDS(kal_adm2, file.path(DIR_DERIVED, "kal_adm2.rds"))
saveRDS(kal_prov, file.path(DIR_DERIVED, "kal_prov.rds"))

# --- (b) FIRMS: archive retention probe -------------------------------------
# Establishes, rather than assumes, where the retrievable record begins. Kept in
# the pipeline so the reported window can be re-verified by a marker.
firms_url <- function(source, bbox, days, start)
  glue("https://firms.modaps.eosdis.nasa.gov/api/area/csv/{map_key}/{source}/",
       "{paste(bbox, collapse=',')}/{days}/{format(start, '%Y-%m-%d')}")

# A single transient failure would otherwise become a silent hole in the daily
# count series, which is indistinguishable from a day with no fire. Retry with
# backoff, and return an attribute recording whether the block was ever served.
firms_get <- function(source, bbox, days, start, tries = 4L) {
  for (k in seq_len(tries)) {
    out <- try(suppressWarnings(read_csv(firms_url(source, bbox, days, start),
                                         show_col_types = FALSE, progress = FALSE)),
               silent = TRUE)
    ok <- !inherits(out, "try-error") && "latitude" %in% names(out)
    if (ok && nrow(out) > 0) return(structure(out, served = TRUE))
    if (ok && nrow(out) == 0 && k == tries) return(structure(tibble(), served = TRUE))
    Sys.sleep(2 * k)
  }
  structure(tibble(), served = FALSE)
}

KBB <- unname(KALIMANTAN_BBOX[c("west","south","east","north")])
probe_dates <- as.Date(c("2026-01-05","2026-03-05","2026-05-05","2026-06-05",
                         "2026-06-30","2026-07-01","2026-07-02","2026-08-05"))
say("probing FIRMS NRT retention over a Kalimantan-wide box")
probe <- map_dfr(probe_dates, function(d) {
  tibble(date = d, n = nrow(firms_get("VIIRS_SNPP_NRT", KBB, 1, d)))
})
# Confirm the SP archive really is empty for 2026 rather than untested.
probe_sp <- nrow(firms_get("VIIRS_SNPP_SP", KBB, 2, as.Date("2026-08-05")))
say("SP archive rows for a 2026-08-05 Kalimantan query: {probe_sp}")
say("retention probe: {paste(sprintf('%s=%d', probe$date, probe$n), collapse='  ')}")
saveRDS(list(nrt = probe, sp_rows = probe_sp), file.path(DIR_DERIVED, "retention_probe.rds"))

# --- (b) FIRMS: the study pull ---------------------------------------------
# NRT only: the SP archive is empty for 2026, verified immediately above.
starts <- seq(WINDOW_START, WINDOW_END, by = paste(FIRMS_MAX_DAY_RANGE, "days"))
pull_platform <- function(source, label) {
  map_dfr(starts, function(s) {
    days <- as.integer(min(WINDOW_END, s + FIRMS_MAX_DAY_RANGE - 1L) - s) + 1L
    d <- firms_get(source, KBB, days, s)
    served <- isTRUE(attr(d, "served"))
    say("  {label} {format(s,'%d %b')} +{days}d: {nrow(d)} rows",
        "{if (!served) '  << NOT SERVED after retries, treat as missing not zero' else ''}")
    blocks <<- add_row(blocks, platform = label, start = s, days = days,
                       n = nrow(d), served = served)
    if (!nrow(d)) return(NULL)
    d |> mutate(platform = label, acq_source = "NRT")
  })
}
blocks <- tibble(platform = character(), start = as.Date(character()),
                 days = integer(), n = integer(), served = logical())
say("pulling FIRMS VIIRS C2 NRT, {WINDOW_START} to {WINDOW_END}, both platforms")
fires_raw <- bind_rows(
  pull_platform("VIIRS_SNPP_NRT",   "S-NPP"),
  pull_platform("VIIRS_NOAA20_NRT", "NOAA-20")
)
say("FIRMS rows downloaded: {format(nrow(fires_raw), big.mark=',')}")
say("columns: {paste(names(fires_raw), collapse=', ')}")
say("observed date range: {min(fires_raw$acq_date)} to {max(fires_raw$acq_date)}")
saveRDS(fires_raw, file.path(DIR_DERIVED, "firms_raw.rds"))
# Block-level provenance: which 5-day requests were served, and which were not.
# A block that was never served is missing data, not an observed zero.
saveRDS(blocks, file.path(DIR_DERIVED, "request_blocks.rds"))
unserved <- blocks |> filter(!served)
if (nrow(unserved)) {
  say("WARNING {nrow(unserved)} request block(s) never served; the daily series ",
      "has genuine gaps at: {paste(unique(format(unserved$start, '%d %b')), collapse=', ')}")
} else say("all {nrow(blocks)} request blocks served")

# --- (b) evidence for the choice of study regency ---------------------------
# Ranked over the whole window, not a convenience sample, and after the same
# confidence screen the analysis will apply.
pts <- fires_raw |>
  filter(!is.na(latitude), !is.na(longitude),
         between(latitude, -90, 90), between(longitude, -180, 180),
         tolower(substr(confidence, 1, 1)) %in% c("n", "h")) |>
  st_as_sf(coords = c("longitude","latitude"), crs = 4326)

area_km2 <- kal_adm2 |> st_transform(STUDY_CRS) |>
  mutate(a = as.numeric(set_units(st_area(geometry), "km^2"))) |>
  st_drop_geometry() |> select(shapeName, province, area_km2 = a)

regency_ranking <- st_join(pts, kal_adm2[, "shapeName"], join = st_within) |>
  st_drop_geometry() |> filter(!is.na(shapeName)) |>
  count(shapeName, name = "n_detections") |>
  left_join(area_km2, by = "shapeName") |>
  mutate(per_1000km2 = n_detections / area_km2 * 1000) |>
  arrange(desc(n_detections))
say("top regency by detections: {regency_ranking$shapeName[1]} ",
    "({format(regency_ranking$n_detections[1], big.mark=',')} detections, ",
    "{round(regency_ranking$per_1000km2[1])} per 1,000 km2)")
say("configured STUDY_REGENCY: {STUDY_REGENCY} ",
    "(rank {which(regency_ranking$shapeName == STUDY_REGENCY)})")
saveRDS(regency_ranking, file.path(DIR_DERIVED, "regency_ranking.rds"))

writeLines(c("# Acquisition log", "",
  glue("Run at: {format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')}"),
  glue("Access date: {ACCESS_DATE}"),
  glue("Window: {WINDOW_START} to {WINDOW_END} ({as.integer(WINDOW_END-WINDOW_START)+1} days)"),
  glue("Sensors: VIIRS 375 m C2 NRT, S-NPP and NOAA-20"),
  glue("sf {packageVersion('sf')}; GDAL {sf_extSoftVersion()[['GDAL']]}; PROJ {sf_extSoftVersion()[['PROJ']]}"),
  "", "## Steps", "", paste0("- ", log_lines)),
  file.path(DIR_DERIVED, "acquisition_log.md"))
say("acquisition complete")
