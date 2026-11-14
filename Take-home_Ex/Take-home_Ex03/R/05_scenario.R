# ---------------------------------------------------------------------------
# 05_scenario.R - the 2030 network: JRL and CRL phase 1 added
#
# The calibrated model is applied twice to the same 10,817 HDB blocks: once to
# today's stations and once to the 2030 set, in which the new stations exist and
# the interchanges they create count their extra line. Comparing the two runs
# of the same model isolates the effect of the new stations; comparing either
# with observed ridership would mix it with the model's own error.
#
# Prediction intervals come from the bootstrap draws of R/04_calibrate.R: each
# draw's parameters are re-scaled to the observed total and re-applied.
# Run from Take-home_Ex/Take-home_Ex03/ after R/04:  Rscript R/05_scenario.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse) })
source("R/00_config.R"); source("R/huff_core.R")

cal <- read_rds(file.path(DER, "calibration.rds"))
stations <- read_rds(file.path(DER, "stations.rds")) |> arrange(name)
si_now <- read_rds(file.path(DER, "si_existing.rds"))
si_30  <- read_rds(file.path(DER, "si_2030.rds"))
best <- cal$best; vars <- cal$models[[best]]; th <- cal$fits[[best]]$theta; k <- cal$fits[[best]]$k

stn_now <- filter(stations, existing)
stn_30  <- stations |> mutate(n_lines = n_lines_2030, lrt_only = lrt_only_2030)

design_for <- function(stn) {
  if (!length(vars)) return(matrix(numeric(0), nrow = nrow(stn), ncol = 0))
  X <- sapply(vars, \(v) if (v == "lrt_only") as.numeric(stn$lrt_only) else log1p(stn[[v]]))
  matrix(X, nrow = nrow(stn), dimnames = list(NULL, vars))
}
pairs_for <- function(si, stn) {
  o_ids <- sort(unique(si$origin_id))
  si |> filter(name %in% stn$name) |>
    transmute(origin_id, name, oi = match(origin_id, o_ids), sj = match(name, stn$name),
              logd = log(dist), du = du)
}
P_now <- pairs_for(si_now, stn_now); X_now <- design_for(stn_now)
P_30  <- pairs_for(si_30, stn_30);   X_30  <- design_for(stn_30)

run <- function(theta, k, P, X, n) k * huff_predict(theta, P, X, n)
pred_now <- run(th, k, P_now, X_now, nrow(stn_now))
pred_30  <- run(th, k, P_30, X_30, nrow(stn_30))

res <- stn_30 |> st_drop_geometry() |>
  select(name, existing, new, lines_now, lines_2030, codes_2030, PLN_AREA_N, REGION_N,
         bus_stops, poi, tap_am) |>
  mutate(pred_2030 = pred_30) |>
  left_join(tibble(name = stn_now$name, pred_now = pred_now), by = "name") |>
  mutate(change = pred_2030 - pred_now, change_pct = change / pred_now,
         line_new = case_when(str_detect(lines_2030, "JRL") & new ~ "JRL",
                              str_detect(lines_2030, "CRL") & new ~ "CRL", TRUE ~ NA))

# Bootstrap intervals ----------------------------------------------------------------
obs <- cal$stn$tap_am
X_cal <- design_for(cal$stn |> arrange(name))
P_cal <- pairs_for(si_now, cal$stn |> arrange(name))
draws <- cal$boot |> pivot_wider(names_from = term, values_from = estimate)
boot_pred <- map(seq_len(nrow(draws)), \(b) {
  th_b <- unlist(draws[b, names(th)])
  a <- huff_predict(th_b, P_cal, X_cal, nrow(X_cal))
  ci_ <- cal$cal_idx
  k_b <- exp(mean(log(obs[ci_]) - log(pmax(a[ci_], 1e-9))))
  tibble(name = stn_30$name, b = b, pred_2030 = run(th_b, k_b, P_30, X_30, nrow(stn_30))) |>
    left_join(tibble(name = stn_now$name, pred_now = run(th_b, k_b, P_now, X_now, nrow(stn_now))), by = "name")
}) |> list_rbind()
pi <- boot_pred |> group_by(name) |>
  summarise(pred_2030_lo = quantile(pred_2030, 0.025), pred_2030_hi = quantile(pred_2030, 0.975),
            change_lo = quantile(pred_2030 - pred_now, 0.025, na.rm = TRUE),
            change_hi = quantile(pred_2030 - pred_now, 0.975, na.rm = TRUE), .groups = "drop")
res <- left_join(res, pi, by = "name")

# Market areas -------------------------------------------------------------------------
# Each block's dominant station (highest share), today and in 2030, and the
# share-weighted dwelling units every station draws.
share_tbl <- function(P, X) P |> mutate(p = huff_shares(th, P, X))
sh_now <- share_tbl(P_now, X_now); sh_30 <- share_tbl(P_30, X_30)
dominant <- function(sh) sh |> group_by(origin_id) |> slice_max(p, n = 1, with_ties = FALSE) |>
  ungroup() |> select(origin_id, station = name, p_max = p)
dom <- dominant(sh_now) |> rename(station_now = station, p_now = p_max) |>
  left_join(dominant(sh_30) |> rename(station_2030 = station, p_2030 = p_max), by = "origin_id") |>
  mutate(switched = station_now != station_2030,
         to_new = station_2030 %in% res$name[res$new])

catch <- full_join(
  sh_now |> group_by(name) |> summarise(du_now = sum(du * p), .groups = "drop"),
  sh_30 |> group_by(name) |> summarise(du_2030 = sum(du * p), .groups = "drop"), by = "name")
res <- left_join(res, catch, by = "name")

write_rds(list(res = res, dom = dom, boot_pred = boot_pred, sh_now = sh_now, sh_30 = sh_30),
          file.path(DER, "scenario.rds"))
res |> filter(new) |> arrange(-pred_2030) |> select(name, line_new, pred_2030, pred_2030_lo, pred_2030_hi) |> print(n = 40)
res |> filter(existing) |> arrange(change_pct) |> select(name, pred_now, pred_2030, change_pct) |> head(15) |> print()
message(sum(dom$switched), " blocks change dominant station; ", sum(dom$to_new), " now dominated by a new station")
