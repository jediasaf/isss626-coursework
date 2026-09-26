# ---------------------------------------------------------------------------
# 08 - Extended r range
#
# The 5 km cap used in R/04 turned out to be too short to answer the question it
# was set up to answer. L_inhom and g_inhom remain outside their envelopes across
# the whole of that range, so 5 km establishes only that dependence extends at
# least that far, not where it ends.
#
# The r range is therefore extended to 20 km. This is not a post-hoc licence: the
# pre-stated rule in R/00_config.R caps r at a quarter of the window's shortest
# dimension, and the window is 86 km wide, so 20 km was always inside the
# admissible range. The 5 km choice was conservative, and it was wrong.
#
# nsim is reduced to 99 because the pair count grows with r^2. The Monte Carlo
# floor becomes 0.01 rather than 0.005, which is stated wherever it matters.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(spatstat.geom); library(spatstat.explore); library(spatstat.random)
})
source("R/00_config.R")
set.seed(SEED)
P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
FO <- readRDS(file.path(DIR_DERIVED, "firstorder.rds"))
X  <- P$ppp_cell
sigma <- FO$sigma_main

RMAX_EXT <- 20000
NSIM_EXT <- 99L
r_ext <- seq(0, RMAX_EXT, length.out = 321)

lam_data <- density(X, sigma = sigma, edge = TRUE, diggle = TRUE, positive = TRUE)

Linhom_reest <- function(Y, ...) {
  lam <- density(Y, sigma = sigma, edge = TRUE, diggle = TRUE, positive = TRUE)
  a <- list(...); a$correction <- "border"
  do.call(Linhom, c(list(Y, lambda = lam), a))
}
# pcfinhom offers no border correction (only isotropic, translation or none), and
# the translation correction over ~7 million pairs per pattern does not complete in
# a workable time for 99 simulated patterns. The work is therefore split:
#
#   significance over 0 to 20 km  -> L_inhom with border correction and envelopes
#   curve shape over 0 to 20 km   -> a single translation-corrected g_inhom, no
#                                    simulations, which costs one pass
#   envelope on g_inhom           -> the 0 to 3 km version already in R/04
#
# This is stated rather than hidden, because it means no significance claim is made
# for g_inhom beyond 3 km; the claim over the extended range rests on L_inhom.
t0 <- Sys.time()
env_L_ext <- envelope(X, Linhom_reest, r = r_ext, nsim = NSIM_EXT,
                      simulate = expression(rpoispp(lam_data)),
                      savefuns = FALSE, verbose = FALSE)
message(sprintf("Linhom envelope to %d km: %.1f min", RMAX_EXT/1000,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

lam_for_obs <- lam_data
t0 <- Sys.time()
g_obs_ext <- pcfinhom(X, lambda = lam_for_obs, r = r_ext,
                      divisor = "d", correction = "translate")
message(sprintf("observed pcfinhom to %d km (no simulations): %.1f min", RMAX_EXT/1000,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Where does the excess actually end? Two readings.
#
# (a) L_inhom against its envelope: the first lag beyond the peak from which the
#     observed curve stays inside the envelope for all larger lags.
L <- as.data.frame(env_L_ext) |> mutate(cobs = obs - r, chi = hi - r) |> filter(r > 200)
inside <- L$cobs <= L$chi
stays_in <- rev(Reduce(function(a, b) a && b, rev(inside), accumulate = TRUE))
idx <- which(stays_in)[1]
r_L_reenter <- if (is.na(idx)) NA_real_ else L$r[idx]

# (b) the observed g_inhom: where it falls to within 10% of its asymptote of 1,
#     a descriptive statement about the correlation range that does not depend on
#     an envelope.
G <- as.data.frame(g_obs_ext) |> filter(r > 200, is.finite(trans))
g_at <- function(target) G$trans[which.min(abs(G$r - target))]
below <- which(G$trans <= 1.1)
r_g_decay <- if (length(below)) G$r[below[1]] else NA_real_

summary_ext <- tibble(
  peak_r_m = G$r[which.max(G$trans)], peak_g = max(G$trans),
  g_1km = g_at(1000), g_2km = g_at(2000), g_5km = g_at(5000),
  g_10km = g_at(10000), g_20km = g_at(20000),
  r_g_within_10pct_of_1 = r_g_decay,
  r_L_reenters_envelope = r_L_reenter)
print(as.data.frame(summary_ext |> mutate(across(everything(), ~ signif(.x, 4)))))

# Cross-check that the extended curve reproduces the R/04 estimate over the lags
# both cover. Both use the translation correction, so they should agree closely;
# a disagreement would indicate that the extended r grid had changed the estimate.
SO <- readRDS(file.path(DIR_DERIVED, "secondorder.rds"))
g_tr <- as.data.frame(SO$env_pcfinhom)
cmp_r <- seq(300, 2900, by = 100)
cmp <- tibble(r = cmp_r,
              r04 = approx(g_tr$r, g_tr$obs, cmp_r)$y,
              r08 = approx(G$r, G$trans, cmp_r)$y) |>
  mutate(rel_diff = (r08 - r04) / r04)
message(sprintf("cross-check over 0.3-2.9 km: median relative difference %.2f%%, max %.2f%%",
                100 * median(abs(cmp$rel_diff), na.rm = TRUE),
                100 * max(abs(cmp$rel_diff), na.rm = TRUE)))

# fv and envelope objects carry the closures that built them, and those closures
# capture the calling environment, which here includes other multi-megabyte caches.
# Storing plain data frames keeps every number and drops hundreds of megabytes.
slim_fv <- function(x) if (inherits(x, c("fv", "envelope"))) as.data.frame(x) else x
saveRDS(lapply(list(r_ext = r_ext, nsim = NSIM_EXT, rmax = RMAX_EXT,
             env_L_ext = env_L_ext, g_obs_ext = g_obs_ext,
             summary_ext = summary_ext, correction_cmp = cmp),
        slim_fv), file.path(DIR_DERIVED, "extended.rds"))
message("extended-range analysis cached")
