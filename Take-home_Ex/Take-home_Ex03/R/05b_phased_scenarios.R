# ---------------------------------------------------------------------------
# 05b_phased_scenarios.R - networks defined by official opening schedules
#
# R/05 runs the comparison the brief asks for: today's network against today's
# plus the JRL and the CRL phase 1. The brief calls that the network "by 2030";
# whether it really is depends on LTA's current schedules, which are recorded,
# with sources, in future_station_schedule.csv. This script builds three
# networks from that table and runs the same calibrated model on each:
#   brief       today + JRL + CRL phase 1 (the brief's comparison, as in R/05)
#   end_2030    today + every planned station LTA currently schedules to open
#               by the end of 2030
#   full        today + every planned station in the station list
# All scenarios hold today's homes and today's surroundings fixed, so they show
# where today's HDB residents would board, not ridership in a given year.
# Run from Take-home_Ex/Take-home_Ex03/ after R/04:  Rscript R/05b_phased_scenarios.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse) })
source("R/00_config.R"); source("R/prep_functions.R"); source("R/huff_core.R")
sf_use_s2(FALSE)

cal <- read_rds(file.path(DER, "calibration.rds"))
origins <- read_rds(file.path(DER, "origins.rds"))
codes_all <- read_rds(file.path(DER, "station_codes.rds"))
subzone <- st_read(F_SUBZONE, quiet = TRUE) |> st_transform(STUDY_CRS) |> st_make_valid() |>
  select(SUBZONE_N, PLN_AREA_N, REGION_N)
poi_sf <- st_read(F_POI, quiet = TRUE) |> st_transform(STUDY_CRS) |> filter(group %in% POI_GROUPS)
busstop <- if (file.exists(F_BUSSTOP)) {
  shp <- grep("\\.shp$", unzip(F_BUSSTOP, list = TRUE)$Name, value = TRUE)[1]
  st_read(paste0("/vsizip/", F_BUSSTOP, "/", shp), quiet = TRUE) |> st_set_crs(STUDY_CRS)
} else st_read(file.path(DER, "osm_busstops_sg.gpkg"), quiet = TRUE) |> st_transform(STUDY_CRS)

sched <- read_csv(F_SCHEDULE, show_col_types = FALSE, comment = "#")
planned <- codes_all |> filter(scenario != "existing") |>
  left_join(select(sched, station_code, stage, planned_opening), by = "station_code")
missing_sched <- planned |> filter(is.na(planned_opening)) |> pull(station_code)
if (length(missing_sched)) message("No verified opening year (kept out of end_2030): ", paste(missing_sched, collapse = ", "))

scenario_codes <- list(
  brief    = c(codes_all$station_code[codes_all$scenario == "existing"],
               planned$station_code[planned$scenario %in% c("JRL", "CRL1")]),
  end_2030 = c(codes_all$station_code[codes_all$scenario == "existing"],
               planned$station_code[!is.na(planned$planned_opening) & planned$planned_opening <= 2030]),
  full     = codes_all$station_code)

existing_names <- unique(codes_all$station_name[codes_all$scenario == "existing"])
vars <- cal$models[[cal$best]]; th <- cal$fits[[cal$best]]$theta; k <- cal$fits[[cal$best]]$k
design <- function(s) {
  X <- sapply(vars, \(v) if (v == "lrt_only") as.numeric(s$lrt_only) else log1p(s[[v]]))
  matrix(X, nrow = nrow(s), dimnames = list(NULL, vars))
}
run_network <- function(code_set) {
  stn <- complex_of(filter(codes_all, station_code %in% code_set)) |> add_attrac() |> arrange(name)
  si <- si_frame(origins, stn)
  o_ids <- sort(unique(si$origin_id))
  P <- si |> transmute(oi = match(origin_id, o_ids), sj = match(name, stn$name), logd = log(dist), du)
  stn |> st_drop_geometry() |>
    transmute(name, lines, n_lines, new = !name %in% existing_names,
              pred = k * huff_predict(th, P, design(stn), nrow(stn)))
}
today <- run_network(codes_all$station_code[codes_all$scenario == "existing"]) |> select(name, pred_today = pred)
runs <- imap_dfr(scenario_codes, \(cs, nm) run_network(cs) |> mutate(scenario = nm)) |>
  left_join(today, by = "name") |>
  mutate(change_pct = pred / pred_today - 1)

summary_tbl <- runs |> group_by(scenario) |>
  summarise(stations = n(), new_stations = sum(new),
            new_boardings = sum(pred[new]), total = sum(pred), .groups = "drop")
write_rds(list(runs = runs, summary = summary_tbl, scenario_codes = scenario_codes,
               planned = planned, sched = sched), file.path(DER, "phased_scenarios.rds"))
print(summary_tbl)
runs |> filter(!new) |> group_by(name) |> summarise(across(change_pct, list(min = min, max = max))) |>
  arrange(change_pct_min) |> head(12) |> print()
