# ---------------------------------------------------------------------------
# 04b_diagnostics.R - checks on the data and the fitted model
#
# Every check here answers a question a reviewer could ask, and every number
# the report quotes about data quality or model validity comes from here:
#   - is any station counted twice, or any passenger-data code left unused?
#   - do the station coordinates agree with LTA's own station polygons?
#   - did the HDB join put two property records on one footprint?
#   - did every optimisation converge, including the bootstrap refits?
#   - is total demand conserved between scenarios (redistribution, not growth)?
#   - are the attraction variables collinear enough to destabilise exponents?
#   - is the Moran's I of the residuals robust to a permutation test?
#   - are the newly opened Circle Line stations still ramping up?
#   - what does the trip rate k mean in everyday units?
# Run from Take-home_Ex/Take-home_Ex03/ after R/04:  Rscript R/04b_diagnostics.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse); library(spdep) })
source("R/00_config.R"); source("R/huff_core.R")
set.seed(SEED)
sf_use_s2(FALSE)

cal      <- read_rds(file.path(DER, "calibration.rds"))
stations <- read_rds(file.path(DER, "stations.rds"))
origins  <- read_rds(file.path(DER, "origins.rds"))
rid      <- read_rds(file.path(DER, "ridership.rds"))
taps_m   <- read_rds(file.path(DER, "taps_by_month.rds"))
D <- list()

# 1. Passenger-data codes: one row per station complex, nothing double counted
D$pv_codes_per_station <- rid |> mutate(n_codes = str_count(pv_codes, fixed(" + ")) + 1) |>
  count(n_codes, name = "stations")
D$pv_codes_total <- n_distinct(taps_m$pt_code)
D$pv_codes_used  <- sum(str_count(rid$pv_codes, fixed(" + ")) + 1)

# 2. Station coordinates against LTA's station polygons (same code)
lta <- st_read(paste0("/vsizip/", F_LTA_STN, "/TrainStation_Sep2026/RapidTransitSystemStation.shp"),
               quiet = TRUE) |> st_set_crs(STUDY_CRS) |> st_make_valid()
codes_pts <- bind_rows(read_csv(F_STN_OP, show_col_types = FALSE),
                       read_csv(F_STN_FUT, show_col_types = FALSE)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |> st_transform(STUDY_CRS)
lta_codes <- lta |> mutate(code = str_split(STN_NAM, "[/ ,]+")) |> unnest(code)
m <- inner_join(st_drop_geometry(codes_pts) |> mutate(i = row_number()),
                st_drop_geometry(lta_codes) |> mutate(k = row_number()), by = c("station_code" = "code"))
m$dist_m <- as.numeric(st_distance(codes_pts[m$i, ], lta_codes[m$k, ], by_element = TRUE))
D$coord_check <- m |> group_by(station_code, station_name) |> summarise(dist_m = min(dist_m), .groups = "drop")
D$coord_summary <- D$coord_check |> summarise(n = n(), median_m = median(dist_m),
                                              p95_m = quantile(dist_m, 0.95), max_m = max(dist_m),
                                              within_100m = mean(dist_m <= 100))

# 3. HDB join: one footprint per property record
D$origins_dup_postal <- origins |> st_drop_geometry() |> filter(placed == "footprint") |>
  count(blk_no, street) |> filter(n > 1) |> nrow()
D$origins_du_total <- sum(origins$du)
D$origins_empty_geom <- sum(st_is_empty(origins))
xy <- st_coordinates(origins)
D$origins_same_xy <- sum(duplicated(round(xy, 1)))

# 4. Convergence of every fit
D$convergence <- tibble(model = names(cal$fits), code = map_int(cal$fits, "convergence"))

# 5. Demand conservation and the meaning of k
scn <- read_rds(file.path(DER, "scenario.rds"))
D$conservation <- scn$res |> summarise(now = sum(pred_now, na.rm = TRUE), y2030 = sum(pred_2030)) |>
  mutate(diff_pct = 100 * (y2030 / now - 1))
k <- cal$fits[[cal$best]]$k
D$k <- k
# k is monthly AM-peak boardings per dwelling unit; per weekday it divides by
# the average number of weekdays in the three months.
# Monday-to-Friday days: 21, 22 and 22. National Day (9 Aug 2026) falls on a
# Sunday, so Monday 10 Aug is the holiday, which DataMall counts with weekends.
wk <- c(`2026-08` = 20, `2026-09` = 22, `2026-10` = 22)
D$k_per_weekday <- k / mean(wk)

# 6. Collinearity of the attraction variables (fitted stations)
stn <- cal$stn[cal$cal_idx, ] |> st_drop_geometry()
X <- tibble(n_lines = log1p(stn$n_lines), bus_stops = log1p(stn$bus_stops),
            poi = log1p(stn$poi), lrt_only = as.numeric(stn$lrt_only))
D$cor_X <- cor(X)
D$vif <- map_dbl(names(X), \(v) 1 / (1 - summary(lm(reformulate(setdiff(names(X), v), v), data = X))$r.squared)) |>
  set_names(names(X))
D$boot_cor <- cal$boot |> pivot_wider(names_from = term, values_from = estimate) |> select(-b) |> cor()

# 7. Moran's I: analytical (randomisation) and permutation
fit_stn <- cal$stn[cal$cal_idx, ]
lw <- nb2listw(knn2nb(knearneigh(st_coordinates(fit_stn), k = 6)), style = "W")
D$moran_perm <- moran.mc(fit_stn$resid_log, lw, nsim = 999)
D$moran_k <- map_dfr(c(4, 6, 8, 10), \(kk) {
  lwk <- nb2listw(knn2nb(knearneigh(st_coordinates(fit_stn), k = kk)), style = "W")
  mt <- moran.mc(fit_stn$resid_log, lwk, nsim = 999)
  tibble(k = kk, I = mt$statistic, p = mt$p.value)
})

# 8. Circle Line stage 6: ramp-up across the three months
D$ccl6 <- taps_m |> filter(str_detect(pt_code, "CC30|CC31|CC32")) |>
  select(YEAR_MONTH, pt_code, tap_am, tap_wkday)

# 9. Residuals against station characteristics (association only)
D$resid_by_group <- fit_stn |> st_drop_geometry() |>
  mutate(group = case_when(lrt_only ~ "LRT only", n_lines > 1 ~ "Interchange",
                           REGION_N == "CENTRAL REGION" ~ "Central, one line", TRUE ~ "Other, one line")) |>
  group_by(group) |>
  summarise(n = n(), median_resid = median(resid_log), mean_resid = mean(resid_log), .groups = "drop")

write_rds(D, file.path(DER, "diagnostics.rds"))
str(D[c("pv_codes_per_station", "pv_codes_total", "pv_codes_used", "coord_summary", "origins_dup_postal",
        "origins_du_total", "origins_same_xy", "convergence", "conservation", "k", "k_per_weekday", "vif")], max.level = 2)
print(round(D$cor_X, 2)); print(round(D$boot_cor, 2)); print(D$moran_perm); print(D$moran_k)
print(D$ccl6); print(D$resid_by_group); print(D$coord_check |> arrange(-dist_m) |> head(8))
