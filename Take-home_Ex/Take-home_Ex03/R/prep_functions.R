# ---------------------------------------------------------------------------
# prep_functions.R - the station and catchment rules, shared by R/03_prepare.R
# and the phased network scenarios (R/05b_phased_scenarios.R), so that every
# scenario is built with exactly the same definitions.
# Expects in the calling environment: origins, subzone, poi_sf, busstop.
# ---------------------------------------------------------------------------

# Which line a station code belongs to. Branches and extensions count as their
# parent line (CG is the East-West Line's Changi branch, CE the Circle Line's
# extension), so a station on both is not counted as an interchange.
line_of <- function(code) {
  p <- str_extract(code, "^[A-Z]+")
  case_when(p %in% c("EW", "CG") ~ "EW", p %in% c("CC", "CE") ~ "CC",
            p == "BP" ~ "BPLRT", p %in% c("SE", "SW", "STC") ~ "SKLRT",
            p %in% c("PE", "PW", "PTC") ~ "PGLRT", p %in% JRL_PREFIX ~ "JRL",
            p %in% c("CR", "CP") ~ "CRL", TRUE ~ p)
}

# One destination per physical station: codes are grouped by station name.
complex_of <- function(df) df |>
  group_by(name = station_name) |>
  summarise(codes = paste(sort(station_code), collapse = "/"),
            lines = paste(sort(unique(line)), collapse = "+"),
            n_lines = n_distinct(line),
            lrt_only = all(str_detect(line, "LRT")),
            lon = mean(lon), lat = mean(lat), .groups = "drop") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |> st_transform(STUDY_CRS)

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

# One row per origin and candidate station: the stations within
# CHOICE_RADIUS_M, topped up to the K_MIN nearest.
si_frame <- function(orig, stn, radius = CHOICE_RADIUS_M) {
  D <- st_distance(orig, stn) |> units::drop_units()
  keep <- D <= radius
  r <- t(apply(D, 1, rank, ties.method = "first")) <= K_MIN
  idx <- which(keep | r, arr.ind = TRUE)
  tibble(origin_id = orig$origin_id[idx[, 1]], name = stn$name[idx[, 2]],
         dist = pmax(D[idx], MIN_DIST_M)) |>
    left_join(st_drop_geometry(select(orig, origin_id, du)), by = "origin_id")
}
