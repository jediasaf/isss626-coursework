# ---------------------------------------------------------------------------
# 06_value_added.R - does distance decay differ across Singapore?
#
# The calibrated Huff model has one distance-decay exponent for the whole
# island. That is the global-model assumption Lessons 9 and 10 questioned for
# regression: one lambda says a resident of Punggol and a resident of Bukit
# Merah trade distance against station quality the same way. Feeder bus
# networks, block layouts and walking environments differ by region, and the
# chosen model's residuals are spatially autocorrelated (R/04), so the decay
# may differ too.
#
# Two extensions of the chosen model, scored by the same 10-fold station
# cross-validation as R/04:
#   R5  one lambda per planning region of the origin block (five regions)
#   EXP the same attraction terms with exponential decay exp(beta * d_km),
#       the usual alternative to the power function
# Run from Take-home_Ex/Take-home_Ex03/ after R/04:  Rscript R/06_value_added.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(tidyverse); library(spdep) })
source("R/00_config.R"); source("R/huff_core.R")
set.seed(SEED)

cal <- read_rds(file.path(DER, "calibration.rds"))
origins <- read_rds(file.path(DER, "origins.rds"))
si <- read_rds(file.path(DER, "si_existing.rds"))
stn <- cal$stn; obs <- stn$tap_am; cal_idx <- cal$cal_idx
vars <- cal$models[[cal$best]]
X <- sapply(vars, \(v) if (v == "lrt_only") as.numeric(stn$lrt_only) else log1p(stn[[v]]))
X <- matrix(X, nrow = nrow(stn), dimnames = list(NULL, vars))

regions <- sort(unique(origins$REGION_N))
pairs <- si |> filter(name %in% stn$name) |>
  left_join(st_drop_geometry(select(origins, origin_id, REGION_N)), by = "origin_id") |>
  transmute(oi = match(origin_id, cal$o_ids), sj = match(name, stn$name),
            logd = log(dist), d_km = dist / 1000, du, reg = factor(REGION_N, regions))

R <- model.matrix(~ reg - 1, pairs)
colnames(R) <- paste0("lambda_", str_to_title(sub(" REGION", "", regions)))
Z_R5  <- cbind(R * pairs$logd, X[pairs$sj, , drop = FALSE])
Z_EXP <- cbind(beta_km = pairs$d_km, X[pairs$sj, , drop = FALSE])
base <- cal$fits[[cal$best]]$theta
specs <- list(R5 = list(Z = Z_R5, start = c(rep(base[["lambda"]], length(regions)), base[-1])),
              EXP = list(Z = Z_EXP, start = c(-1, base[-1])))

fits <- map(specs, \(s) fit_Z(s$Z, pairs, obs, cal_idx, s$start))
insample <- imap_dfr(fits, \(f, nm) fit_stats(obs[cal_idx], f$pred[cal_idx]) |> mutate(model = nm))
cv_pred <- imap(specs, \(s, nm) {
  pred <- rep(NA_real_, length(obs))
  for (f in 1:10) {
    tr <- which(cal$folds != f); te <- which(cal$folds == f)
    pred[te] <- fit_Z(s$Z, pairs, obs, tr, fits[[nm]]$theta)$pred[te]
  }
  pred
})
cv <- imap_dfr(cv_pred, \(p, nm) fit_stats(obs[cal_idx], p[cal_idx]) |> mutate(model = nm))

# Bootstrap the regional lambdas (stations resampled), to see whether the
# differences exceed their uncertainty.
boot_r5 <- map_dfr(1:100, \(b) {
  f <- fit_Z(Z_R5, pairs, obs, sample(cal_idx, replace = TRUE), fits$R5$theta)
  tibble(b = b, term = names(f$theta), estimate = f$theta)
})
lambda_region <- tibble(term = colnames(R), region = regions,
                        lambda = fits$R5$theta[colnames(R)]) |>
  left_join(boot_r5 |> group_by(term) |>
              summarise(lo = quantile(estimate, 0.025), hi = quantile(estimate, 0.975)), by = "term")

# Does the regional model remove the residual clustering?
res5 <- log(obs[cal_idx]) - log(fits$R5$pred[cal_idx])
nb <- knn2nb(knearneigh(st_coordinates(stn[cal_idx, ]), k = 6))
moran_r5 <- moran.test(res5, nb2listw(nb, style = "W"))

write_rds(list(fits = fits, cv = cv, insample = insample, lambda_region = lambda_region,
               boot_r5 = boot_r5, moran_r5 = moran_r5, cv_pred = cv_pred),
          file.path(DER, "value_added.rds"))
print(insample); print(cv); print(lambda_region); print(round(fits$EXP$theta, 3))
message("Moran's I of R5 log residuals: ", round(moran_r5$estimate[1], 3), " (p = ", signif(moran_r5$p.value, 3), ")")
