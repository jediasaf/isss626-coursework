# ---------------------------------------------------------------------------
# 08_sensitivity.R - how much do the results depend on my choices?
#
# 1. One comparison table for every model (M0 to M3, the regional R5 and the
#    exponential EXP), all scored on the same cross-validation folds, with
#    errors in boardings as well as on the log scale.
# 2. The exclusion rule. Stations with fewer than CAL_MIN_DU_1KM HDB dwelling
#    units within 1 km are left out of the fit. Refit with no threshold and with
#    500 and 2,000, and see what changes in the parameters, the validation and
#    the predictions for the new stations.
# 3. The choice set. Refit with a 2 km and a 5 km radius instead of 3 km.
# Run from Take-home_Ex/Take-home_Ex03/ after R/07:  Rscript R/08_sensitivity.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse) })
source("R/00_config.R"); source("R/huff_core.R")
set.seed(SEED)

cal <- read_rds(file.path(DER, "calibration.rds"))
va  <- read_rds(file.path(DER, "value_added.rds"))
origins  <- read_rds(file.path(DER, "origins.rds"))
stations <- read_rds(file.path(DER, "stations.rds")) |> arrange(name)
stn <- cal$stn; obs <- stn$tap_am; vars <- cal$models[[cal$best]]

# Error measures in boardings as well as on the log scale ---------------------------
metrics <- function(o, p) tibble(
  R2_log    = 1 - sum((log(o) - log(p))^2) / sum((log(o) - mean(log(o)))^2),
  MAE       = mean(abs(p - o)),
  RMSE      = sqrt(mean((p - o)^2)),
  WAPE_pct  = 100 * sum(abs(p - o)) / sum(o),
  bias_pct  = 100 * (sum(p) / sum(o) - 1),
  MdAPE_pct = 100 * median(abs(p - o) / o),
  spearman  = cor(o, p, method = "spearman"))

ci <- cal$cal_idx
comparison <- bind_rows(
  metrics(obs[ci], cal$cv_m0[ci]) |> mutate(model = "M0"),
  imap_dfr(cal$cv_pred, \(p, nm) metrics(obs[ci], p[ci]) |> mutate(model = nm)),
  imap_dfr(va$cv_pred, \(p, nm) metrics(obs[ci], p[ci]) |> mutate(model = nm))) |>
  relocate(model)

# Shared machinery for refits ---------------------------------------------------------
design <- function(s) {
  X <- sapply(vars, \(v) if (v == "lrt_only") as.numeric(s$lrt_only) else log1p(s[[v]]))
  matrix(X, nrow = nrow(s), dimnames = list(NULL, vars))
}
si_frame <- function(orig, stn_sf, radius) {
  D <- st_distance(orig, stn_sf) |> units::drop_units()
  r <- t(apply(D, 1, rank, ties.method = "first")) <= K_MIN
  idx <- which(D <= radius | r, arr.ind = TRUE)
  tibble(origin_id = orig$origin_id[idx[, 1]], name = stn_sf$name[idx[, 2]],
         dist = pmax(D[idx], MIN_DIST_M), du = orig$du[idx[, 1]])
}
pairs_for <- function(si, s) {
  o_ids <- sort(unique(si$origin_id))
  si |> filter(name %in% s$name) |>
    transmute(oi = match(origin_id, o_ids), sj = match(name, s$name), logd = log(dist), du)
}
cv_fit <- function(P, X, o, fit_idx) {
  folds <- rep(NA_integer_, length(o)); folds[fit_idx] <- sample(rep(1:10, length.out = length(fit_idx)))
  pred <- rep(NA_real_, length(o))
  for (f in 1:10) {
    tr <- which(!is.na(folds) & folds != f)
    te <- which(folds == f)
    pred[te] <- huff_fit(P, X, o, fit_idx = tr)$pred[te]
  }
  pred
}
stn_now <- filter(stations, existing, !is.na(tap_am), tap_am > 0)
stn_30  <- stations |> mutate(n_lines = n_lines_2030, lrt_only = lrt_only_2030)
new_names <- stations$name[stations$new]
predict_new <- function(theta, k, si30) {
  P <- pairs_for(si30, stn_30)
  tibble(name = stn_30$name, pred = k * huff_predict(theta, P, design(stn_30), nrow(stn_30))) |>
    filter(name %in% new_names)
}
si_now_3k <- read_rds(file.path(DER, "si_existing.rds"))
si_30_3k  <- read_rds(file.path(DER, "si_2030.rds"))
base_new <- predict_new(cal$fits[[cal$best]]$theta, cal$fits[[cal$best]]$k, si_30_3k)

run_variant <- function(label, si_now, si_30, min_du) {
  P <- pairs_for(si_now, stn_now); X <- design(stn_now); o <- stn_now$tap_am
  alloc0 <- huff_predict(c(-2, rep(0, ncol(X))), P, X, nrow(X))
  idx <- which(stn_now$du_1km >= min_du & alloc0 > 0)
  f <- huff_fit(P, X, o, fit_idx = idx, start = cal$fits[[cal$best]]$theta)
  cvp <- cv_fit(P, X, o, idx)
  nw <- predict_new(f$theta, f$k, si_30) |> left_join(base_new, by = "name", suffix = c("", "_base"))
  # Accuracy on the 148 stations of the main fit, so variants are compared on one set
  common <- which(stn_now$name %in% stn$name[ci])
  list(label = label, n_fit = length(idx), theta = f$theta,
       cv_own = metrics(o[idx], cvp[idx]),
       cv_common = metrics(o[intersect(idx, common)], cvp[intersect(idx, common)]),
       new_total_ratio = sum(nw$pred) / sum(nw$pred_base),
       new_rank_rho = cor(nw$pred, nw$pred_base, method = "spearman"),
       new = nw)
}

# 2. Exclusion thresholds --------------------------------------------------------------------
excl <- map(c(0, 500, 1000, 2000), \(m) run_variant(paste0("min ", m, " DU within 1 km"), si_now_3k, si_30_3k, m))

# How far off are the excluded stations under the main model?
excluded <- stn |> st_drop_geometry() |> filter(!in_fit) |>
  transmute(name, PLN_AREA_N, du_1km, observed = tap_am, predicted = pred, ratio = pred / tap_am)

# 3. Choice-set radius -----------------------------------------------------------------------------
radius <- map(c(2000, 5000), \(r) {
  run_variant(paste0(r / 1000, " km choice set"),
              si_frame(origins, stn_now, r), si_frame(origins, stn_30, r), CAL_MIN_DU_1KM)
})

summarise_variant <- function(v) tibble(
  variant = v$label, n_fit = v$n_fit,
  lambda = v$theta[["lambda"]], n_lines = v$theta[["n_lines"]], poi = v$theta[["poi"]],
  lrt_only = v$theta[["lrt_only"]], bus_stops = v$theta[["bus_stops"]],
  cv_R2_log = v$cv_own$R2_log, cv_WAPE_pct = v$cv_own$WAPE_pct, cv_MdAPE_pct = v$cv_own$MdAPE_pct,
  common_R2_log = v$cv_common$R2_log,
  new_total_vs_main = v$new_total_ratio, new_rank_rho = v$new_rank_rho)
sens <- bind_rows(map_dfr(excl, summarise_variant), map_dfr(radius, summarise_variant))

write_rds(list(comparison = comparison, sens = sens, excluded = excluded,
               excl = excl, radius = radius), file.path(DER, "sensitivity.rds"))
print(comparison, width = 200); print(sens, width = 200)
print(excluded |> arrange(ratio), n = 40)
