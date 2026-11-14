# ---------------------------------------------------------------------------
# 03_prepare.R - origins, destinations and the spatial interaction data frames
#
# Builds, in data/derived/:
#   origins.rds        HDB residential blocks with dwelling units (local market potential)
#   stations.rds       station complexes, existing and planned, with attraction variables
#   ridership.rds      observed weekday AM-peak and all-day tap-ins per existing station
#   si_existing.rds    origin x destination pairs for today's network
#   si_2030.rds        the same with the JRL and CRL phase 1 stations added
#
# Run from Take-home_Ex/Take-home_Ex03/ after 00 to 02:  Rscript R/03_prepare.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse) })
source("R/00_config.R")
sf_use_s2(FALSE)

subzone <- st_read(F_SUBZONE, quiet = TRUE) |> st_transform(STUDY_CRS) |> st_make_valid() |>
  select(SUBZONE_N, PLN_AREA_N, REGION_N)

# 1. Origins -----------------------------------------------------------------------
bldg <- st_read(F_BLDG, quiet = TRUE) |> st_transform(STUDY_CRS) |> st_make_valid()
prop <- read_csv(F_PROP, show_col_types = FALSE)
codes <- read_csv(file.path(DER, "street_codes.csv"), show_col_types = FALSE, col_types = cols(.default = "c"))
pts <- read_csv(file.path(DER, "block_points.csv"), show_col_types = FALSE, col_types = cols(.default = "c"))

res <- prop |> filter(residential == "Y", total_dwelling_units >= MIN_DU) |>
  left_join(select(codes, street, ST_COD), by = "street")
on_fp <- res |> inner_join(select(bldg, BLK_NO, ST_COD, POSTAL_COD, geometry),
                           by = c("blk_no" = "BLK_NO", "ST_COD")) |>
  st_as_sf() |> st_point_on_surface() |> mutate(placed = "footprint")
as_pt <- res |> anti_join(st_drop_geometry(on_fp), by = c("blk_no", "street")) |>
  inner_join(pts, by = c("blk_no", "street")) |>
  mutate(POSTAL_COD = postal, x = as.numeric(x), y = as.numeric(y)) |>
  st_as_sf(coords = c("x", "y"), crs = STUDY_CRS) |> mutate(placed = "geocoded point")

origins <- bind_rows(on_fp, as_pt) |>
  transmute(origin_id = paste(blk_no, street, sep = " | "), blk_no, street,
            town = bldg_contract_town, year_completed, du = total_dwelling_units, placed) |>
  st_join(subzone, join = st_nearest_feature)
stopifnot(!anyDuplicated(origins$origin_id))
write_rds(origins, file.path(DER, "origins.rds"))
message(nrow(origins), " origins, ", sum(origins$du), " dwelling units")

# 2. Observed ridership ---------------------------------------------------------------
# LTA's Passenger Volume by Train Stations: tap-ins and tap-outs per station
# code, day type and hour, one file per month. A station's tap-ins are its
# boardings, which is what the model predicts.
pv_files <- Sys.glob(file.path(RAW, PV_GLOB))
if (length(pv_files) == 0) stop("No station passenger volume file in data/raw: run R/00_download.R with DATAMALL_KEY set")
pv <- map_dfr(pv_files, \(z) read_csv(z, show_col_types = FALSE,
              col_types = cols(.default = "c", TOTAL_TAP_IN_VOLUME = "d", TOTAL_TAP_OUT_VOLUME = "d")))
# All months DataMall serves are kept and averaged: one month can carry a
# disruption or a school holiday, and the month-to-month agreement is itself a
# data-quality check (reported in the technical report).
od_month <- paste(range(pv$YEAR_MONTH), collapse = " to ")
day_keep <- DAY_TYPE
taps_month <- pv |>
  mutate(hour = as.integer(TIME_PER_HOUR), wk = DAY_TYPE == day_keep) |>
  group_by(YEAR_MONTH, pt_code = PT_CODE) |>
  summarise(tap_am = sum(TOTAL_TAP_IN_VOLUME[wk & hour %in% AM_HOURS]),
            tapout_am = sum(TOTAL_TAP_OUT_VOLUME[wk & hour %in% AM_HOURS]),
            tap_wkday = sum(TOTAL_TAP_IN_VOLUME[wk]),
            tap_all = sum(TOTAL_TAP_IN_VOLUME), .groups = "drop")
write_rds(taps_month, file.path(DER, "taps_by_month.rds"))
taps <- taps_month |> group_by(pt_code) |>
  summarise(across(c(starts_with("tap_"), tapout_am), \(x) sum(x) / n_distinct(pv$YEAR_MONTH)),
            .groups = "drop")

# 3. Station complexes -------------------------------------------------------------
# One destination per physical station: an interchange is one place to a rider
# walking from home, whatever the number of lines. Codes are grouped by name.
line_of <- function(code) {
  p <- str_extract(code, "^[A-Z]+")
  case_when(p %in% c("EW", "CG") ~ "EW", p %in% c("CC", "CE") ~ "CC",
            p == "BP" ~ "BPLRT", p %in% c("SE", "SW", "STC") ~ "SKLRT",
            p %in% c("PE", "PW", "PTC") ~ "PGLRT", p %in% JRL_PREFIX ~ "JRL",
            p %in% c("CR", "CP") ~ "CRL", TRUE ~ p)
}
op  <- read_csv(F_STN_OP, show_col_types = FALSE) |> mutate(status = "operational")
fut <- read_csv(F_STN_FUT, show_col_types = FALSE) |> mutate(status = "planned")
od_codes <- unique(unlist(strsplit(taps$pt_code, "/")))
# Planned stations that the passenger data shows are already running are moved
# to the existing network (e.g. if Circle Line stage 6 has opened).
now_open <- intersect(setdiff(od_codes, op$station_code), fut$station_code)
# Codes in the passenger data that neither list has (Circle Line stage 6
# platforms at Marina Bay and Bayfront) are attached to the station their
# joint code already names, so the interchange counts the extra line.
known <- c(op$station_code, fut$station_code)
extra <- map_dfr(setdiff(od_codes, known), \(cd) {
  grp <- strsplit(taps$pt_code[str_detect(taps$pt_code, paste0("(^|/)", cd, "($|/)"))][1], "/")[[1]]
  host <- op |> filter(station_code %in% grp) |> slice(1)
  if (nrow(host)) host |> mutate(station_code = cd) else NULL
})
codes_all <- bind_rows(op, extra, filter(fut, station_code %in% now_open) |> mutate(status = "operational"),
                       filter(fut, !station_code %in% now_open) |> mutate(status = "planned")) |>
  mutate(line = line_of(station_code),
         scenario = case_when(status == "operational" ~ "existing",
                              str_extract(station_code, "^[A-Z]+") %in% JRL_PREFIX ~ "JRL",
                              station_code %in% CRL1_CODES ~ "CRL1", TRUE ~ "later"))

complex_of <- function(df) df |>
  group_by(name = station_name) |>
  summarise(codes = paste(sort(station_code), collapse = "/"),
            lines = paste(sort(unique(line)), collapse = "+"),
            n_lines = n_distinct(line),
            lrt_only = all(str_detect(line, "LRT")),
            lon = mean(lon), lat = mean(lat), .groups = "drop") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |> st_transform(STUDY_CRS)

stn_now  <- complex_of(filter(codes_all, scenario == "existing")) |> mutate(new = FALSE)
stn_2030 <- complex_of(filter(codes_all, scenario %in% c("existing", "JRL", "CRL1"))) |>
  mutate(new = !name %in% stn_now$name)
unmatched <- setdiff(od_codes, codes_all$station_code[codes_all$scenario == "existing"])
if (length(unmatched)) warning("Codes in passenger data with no station: ", paste(unmatched, collapse = ", "))

# Observed tap-ins per complex: the passenger data codes interchanges jointly
# ("NS1/EW24"), so each data code is mapped to the complex of the first of
# its parts that the station list knows.
code2name <- setNames(codes_all$station_name, codes_all$station_code)
first_known <- function(code) {
  nm <- code2name[strsplit(code, "/")[[1]]]
  nm[!is.na(nm)][1] %||% NA_character_
}
ridership <- taps |> mutate(name = map_chr(pt_code, first_known)) |>
  filter(!is.na(name)) |> group_by(name) |>
  summarise(across(c(starts_with("tap_"), tapout_am), sum), pv_codes = paste(pt_code, collapse = " + "),
            .groups = "drop") |>
  mutate(month = od_month)
write_rds(ridership, file.path(DER, "ridership.rds"))

# 4. Attraction variables -----------------------------------------------------------
poi_sf <- st_read(F_POI, quiet = TRUE) |> st_transform(STUDY_CRS) |> filter(group %in% POI_GROUPS)
busstop <- if (file.exists(F_BUSSTOP)) {
  shp <- grep("\\.shp$", unzip(F_BUSSTOP, list = TRUE)$Name, value = TRUE)[1]
  st_read(paste0("/vsizip/", F_BUSSTOP, "/", shp), quiet = TRUE) |> st_set_crs(STUDY_CRS)
} else {
  warning("LTA bus stop file missing: using OpenStreetMap bus stops")
  st_read(file.path(DER, "osm_busstops_sg.gpkg"), quiet = TRUE) |> st_transform(STUDY_CRS)
}
count_within <- function(stn, x, d) lengths(st_is_within_distance(stn, x, d))
du_within <- function(stn, d) sapply(st_is_within_distance(stn, origins, d), \(i) sum(origins$du[i]))
add_attrac <- function(stn) {
  stn |> mutate(
    du_1km = du_within(stn, 1000),
    bus_stops = count_within(stn, busstop, CATCH_M),
    poi = count_within(stn, poi_sf, CATCH_M),
    poi_office = count_within(stn, filter(poi_sf, group == "office"), CATCH_M),
    poi_retail_food = count_within(stn, filter(poi_sf, group %in% c("retail", "food")), CATCH_M)) |>
    st_join(subzone, join = st_nearest_feature)
}
# Line counts differ between today and 2030 at interchanges the new lines
# reach (Jurong East, Choa Chu Kang, Boon Lay, Pasir Ris, Hougang, Ang Mo Kio,
# Bright Hill): the model is calibrated on today's values and predicts with 2030's.
now_lines <- st_drop_geometry(stn_now) |>
  select(name, codes_now = codes, lines_now = lines, n_lines_now = n_lines, lrt_only_now = lrt_only)
stations <- add_attrac(stn_2030) |>
  rename(codes_2030 = codes, lines_2030 = lines, n_lines_2030 = n_lines, lrt_only_2030 = lrt_only) |>
  left_join(now_lines, by = "name") |>
  mutate(existing = name %in% stn_now$name,
         n_lines = n_lines_now, lrt_only = lrt_only_now) |>
  left_join(ridership, by = "name")
write_rds(stations, file.path(DER, "stations.rds"))

# 5. Spatial interaction data frames -------------------------------------------------
# One row per origin and candidate station: the stations within CHOICE_RADIUS_M,
# topped up to the K_MIN nearest.
si_frame <- function(orig, stn) {
  D <- st_distance(orig, stn) |> units::drop_units()
  keep <- D <= CHOICE_RADIUS_M
  r <- t(apply(D, 1, rank, ties.method = "first")) <= K_MIN
  idx <- which(keep | r, arr.ind = TRUE)
  tibble(origin_id = orig$origin_id[idx[, 1]], name = stn$name[idx[, 2]],
         dist = pmax(D[idx], MIN_DIST_M)) |>
    left_join(st_drop_geometry(select(orig, origin_id, du)), by = "origin_id")
}
write_rds(si_frame(origins, filter(stations, existing)), file.path(DER, "si_existing.rds"))
write_rds(si_frame(origins, stations), file.path(DER, "si_2030.rds"))
message("done: ", sum(stations$existing), " existing and ", sum(stations$new), " new station complexes")
