# ---------------------------------------------------------------------------
# 04_calibrate.R - calibrate the Huff model on observed AM-peak boardings
#
# Four nested models, each scored in-sample and by 10-fold cross-validation
# over stations. A held-out station stays in every choice set; only its
# ridership is hidden from the fit. That is the situation of a new station,
# whose attributes are known and whose ridership is not.
#   M0  textbook default: every station equally attractive, lambda = -2 (not fitted)
#   M1  distance only: lambda fitted
#   M2  + lines served
#   M3  + bus stops, points of interest and an LRT-only indicator
# Only stations with an HDB catchment (CAL_MIN_DU_1KM) are scored; all
# stations compete in every choice set. The chosen model is checked against
# MCI::huff.shares(), bootstrapped over stations for parameter uncertainty, and
# its residuals tested for spatial autocorrelation.
#
# Run from Take-home_Ex/Take-home_Ex03/ after R/03_prepare.R:  Rscript R/04_calibrate.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse); library(spdep) })
source("R/00_config.R"); source("R/huff_core.R")
set.seed(SEED)

stations <- read_rds(file.path(DER, "stations.rds"))
si <- read_rds(file.path(DER, "si_existing.rds"))

stn <- stations |> filter(existing, !is.na(tap_am), tap_am > 0) |> arrange(name)
o_ids <- sort(unique(si$origin_id))
pairs <- si |> filter(name %in% stn$name) |>
  transmute(oi = match(origin_id, o_ids), sj = match(name, stn$name), logd = log(dist), du = du)
obs <- stn$tap_am

# Stations the residential model is judged on (see CAL_MIN_DU_1KM).
cal_idx <- which(stn$du_1km >= CAL_MIN_DU_1KM)
stn$in_fit <- seq_len(nrow(stn)) %in% cal_idx
message(length(cal_idx), " of ", nrow(stn), " stations enter the fit")

design <- function(vars, s = stn) {
  if (!length(vars)) return(matrix(numeric(0), nrow = nrow(s), ncol = 0))
  X <- sapply(vars, \(v) if (v == "lrt_only") as.numeric(s$lrt_only) else log1p(s[[v]]))
  matrix(X, nrow = nrow(s), dimnames = list(NULL, vars))
}
MODELS <- list(M1 = character(), M2 = "n_lines",
               M3 = c("n_lines", "bus_stops", "poi", "lrt_only"))
stats_on <- function(pred, idx = cal_idx) fit_stats(obs[idx], pred[idx])

# In-sample fits ---------------------------------------------------------------------
fits <- map(MODELS, \(v) huff_fit(pairs, design(v), obs, fit_idx = cal_idx))
m0_alloc <- huff_predict(-2, pairs, design(character()), length(obs))
m0_pred <- exp(mean(log(obs[cal_idx]) - log(m0_alloc[cal_idx]))) * m0_alloc
insample <- bind_rows(stats_on(m0_pred) |> mutate(model = "M0"),
                      imap_dfr(fits, \(f, nm) stats_on(f$pred) |> mutate(model = nm)))

# Cross-validation ---------------------------------------------------------------------
folds <- rep(NA_integer_, nrow(stn))
folds[cal_idx] <- sample(rep(1:10, length.out = length(cal_idx)))
cv_pred <- map(MODELS, \(v) {
  X <- design(v); pred <- rep(NA_real_, length(obs))
  for (f in 1:10) {
    tr <- which(folds != f); te <- which(folds == f)
    pred[te] <- huff_fit(pairs, X, obs, fit_idx = tr)$pred[te]
  }
  pred
})
cv_m0 <- rep(NA_real_, length(obs))
for (f in 1:10) {
  tr <- which(folds != f); te <- which(folds == f)
  cv_m0[te] <- exp(mean(log(obs[tr]) - log(m0_alloc[tr]))) * m0_alloc[te]
}
cv <- bind_rows(stats_on(cv_m0) |> mutate(model = "M0"),
                imap_dfr(cv_pred, \(p, nm) stats_on(p) |> mutate(model = nm)))

params <- imap_dfr(fits, \(f, nm) tibble(model = nm, term = names(f$theta), estimate = f$theta)) |>
  bind_rows(tibble(model = "M0", term = "lambda", estimate = -2))

# Bootstrap the chosen model --------------------------------------------------------------
best <- "M3"
Xb <- design(MODELS[[best]])
boot <- map_dfr(1:200, \(b) {
  idx <- sample(cal_idx, replace = TRUE)
  f <- huff_fit(pairs, Xb, obs, fit_idx = idx, start = fits[[best]]$theta)
  tibble(b = b, term = names(f$theta), estimate = f$theta)
})
ci <- boot |> group_by(term) |>
  summarise(lo = quantile(estimate, 0.025), hi = quantile(estimate, 0.975), .groups = "drop")

# Check against the MCI package -------------------------------------------------------------
th <- fits[[best]]$theta
stn$A <- exp(as.vector(Xb %*% th[-1]))
one_origin <- si |> filter(origin_id == o_ids[1], name %in% stn$name) |>
  left_join(st_drop_geometry(select(stn, name, A)), by = "name") |> as.data.frame()
assign("one_origin", one_origin, envir = globalenv())  # MCI looks its data up by name
mci <- MCI::huff.shares(one_origin, "origin_id", "name", "A", "dist", gamma = 1,
                        lambda = th[["lambda"]], dtype = "pow")
mine <- tibble(name = stn$name[pairs$sj[pairs$oi == 1]],
               p = huff_shares(th, pairs, Xb)[pairs$oi == 1])
mci_check <- max(abs(mine$p - mci$p_ij[match(mine$name, mci$name)]))

# Residuals: are the misses spatially clustered? ----------------------------------------------
stn <- stn |> mutate(pred = fits[[best]]$pred, cv_pred = cv_pred[[best]],
                     resid_log = log(obs) - log(pred))
fit_stn <- stn[cal_idx, ]
nb <- knn2nb(knearneigh(st_coordinates(fit_stn), k = 6))
moran_resid <- moran.test(fit_stn$resid_log, nb2listw(nb, style = "W"))

write_rds(list(stn = stn, fits = fits, params = params, ci = ci, boot = boot,
               insample = insample, cv = cv, cv_pred = cv_pred, cv_m0 = cv_m0, folds = folds,
               m0_pred = m0_pred, cal_idx = cal_idx, best = best, models = MODELS,
               mci_check = mci_check, moran_resid = moran_resid, o_ids = o_ids),
          file.path(DER, "calibration.rds"))
print(insample); print(cv); print(params, n = 20); print(ci)
message("MCI check, max |difference| in shares: ", signif(mci_check, 3))
message("Moran's I of log residuals: ", round(moran_resid$estimate[1], 3),
        " (p = ", signif(moran_resid$p.value, 3), ")")
