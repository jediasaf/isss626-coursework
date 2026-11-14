# ---------------------------------------------------------------------------
# 07_scenario_r5.R - do the 2030 conclusions survive a regional distance decay?
#
# R/06 found that the distance exponent differs by region, and that the
# regional model (R5) predicts held-out stations better than the single-lambda
# model. The JRL runs through the West, where R5's decay is far steeper than
# the island-wide value, so the 2030 scenario of R/05 is re-run with R5 and the
# two sets of predictions compared.
# Run from Take-home_Ex/Take-home_Ex03/ after R/06:  Rscript R/07_scenario_r5.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse) })
source("R/00_config.R"); source("R/huff_core.R")

cal <- read_rds(file.path(DER, "calibration.rds"))
va  <- read_rds(file.path(DER, "value_added.rds"))
scn <- read_rds(file.path(DER, "scenario.rds"))
origins  <- read_rds(file.path(DER, "origins.rds"))
stations <- read_rds(file.path(DER, "stations.rds")) |> arrange(name)
regions  <- sort(unique(origins$REGION_N))
vars <- cal$models[[cal$best]]
th5 <- va$fits$R5$theta; k5 <- va$fits$R5$k

z_for <- function(si, stn) {
  o_ids <- sort(unique(si$origin_id))
  P <- si |> filter(name %in% stn$name) |>
    left_join(st_drop_geometry(select(origins, origin_id, REGION_N)), by = "origin_id") |>
    transmute(origin_id, oi = match(origin_id, o_ids), sj = match(name, stn$name),
              logd = log(dist), du, reg = factor(REGION_N, regions))
  X <- sapply(vars, \(v) if (v == "lrt_only") as.numeric(stn$lrt_only) else log1p(stn[[v]]))
  X <- matrix(X, nrow = nrow(stn), dimnames = list(NULL, vars))
  R <- model.matrix(~ reg - 1, P)
  list(P = P, Z = cbind(R * P$logd, X[P$sj, , drop = FALSE]), n = nrow(stn))
}
stn_now <- filter(stations, existing)
stn_30  <- stations |> mutate(n_lines = n_lines_2030, lrt_only = lrt_only_2030)
A <- z_for(read_rds(file.path(DER, "si_existing.rds")), stn_now)
B <- z_for(read_rds(file.path(DER, "si_2030.rds")), stn_30)

r5 <- tibble(name = stn_30$name, pred_2030_r5 = k5 * alloc_Z(th5, B$Z, B$P, B$n)) |>
  left_join(tibble(name = stn_now$name, pred_now_r5 = k5 * alloc_Z(th5, A$Z, A$P, A$n)), by = "name") |>
  mutate(change_pct_r5 = pred_2030_r5 / pred_now_r5 - 1)
cmp <- scn$res |> select(name, new, existing, line_new, pred_2030, change_pct) |> left_join(r5, by = "name")

summ <- list(
  rho_new = cor(cmp$pred_2030[cmp$new], cmp$pred_2030_r5[cmp$new], method = "spearman"),
  rho_change = with(filter(cmp, existing), cor(change_pct, change_pct_r5, method = "spearman", use = "complete.obs")),
  ratio_new_jrl = with(filter(cmp, new, line_new == "JRL"), sum(pred_2030_r5) / sum(pred_2030)),
  ratio_new_crl = with(filter(cmp, new, line_new == "CRL"), sum(pred_2030_r5) / sum(pred_2030)))
write_rds(list(cmp = cmp, summ = summ), file.path(DER, "scenario_r5.rds"))
print(summ)
cmp |> filter(new) |> arrange(-pred_2030_r5) |> select(name, line_new, pred_2030, pred_2030_r5) |> print(n = 30)
cmp |> filter(existing) |> arrange(change_pct_r5) |> select(name, change_pct, change_pct_r5) |> head(8) |> print()
