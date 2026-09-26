# ---------------------------------------------------------------------------
# 04 - Second-order properties: clustering, inhibition, or randomness
#
# The central methodological argument of this exercise sits here.
#
# A VIIRS hotspot is not a fire. It is a 375 m pixel in which the sensor
# detected a thermal anomaly during one overpass. One large peat fire therefore
# produces many detections, and the resulting pattern is strongly clustered by
# construction. Testing it against complete spatial randomness is close to
# vacuous: CSR is rejected before any analysis is done, and rejecting it says
# nothing about fire behaviour.
#
# Two things are therefore done instead.
#
#  1. The null model is changed. An inhomogeneous Poisson process with
#     intensity lambda-hat(u) from KDE asks a question worth asking: once the
#     large-scale variation in where fires occur is accounted for, is there
#     residual interaction between events? Kinhom, Linhom and pcfinhom answer
#     that; Kest against CSR does not.
#
#  2. The reference intensity is re-estimated for every simulated pattern
#     rather than held fixed at the data estimate. Holding it fixed lets the
#     estimate absorb the very clustering under test, which biases the envelope
#     outward and makes the test conservative in an uncontrolled way.
#
# The homogeneous statistics are still computed and reported, because the
# contrast between the two nulls is itself the finding.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(spatstat.geom); library(spatstat.explore); library(spatstat.random)
})
source("R/00_config.R")
set.seed(SEED)
P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
FO <- readRDS(file.path(DIR_DERIVED, "firstorder.rds"))

X  <- P$ppp_cell
XA <- P$ppp_all
sigma <- FO$sigma_main

# --- r ranges ---------------------------------------------------------------
# Capped well below the window's shorter dimension (86 km). Beyond roughly a
# quarter of that, edge-corrected estimates rest on very few pairs and the
# window's elongated shape dominates the statistic.
RMAX_K   <- 5000    # m, K / L
RMAX_GF  <- 2000    # m, nearest-neighbour and empty-space
RMAX_PCF <- 3000    # m, pair correlation
r_K   <- seq(0, RMAX_K,   length.out = 257)
r_GF  <- seq(0, RMAX_GF,  length.out = 257)
r_pcf <- seq(0, RMAX_PCF, length.out = 257)

NSIM <- as.integer(Sys.getenv("NSIM_2ND", NSIM_ENVELOPE))
message(sprintf("n = %d, nsim = %d", npoints(X), NSIM))

# --- 1. homogeneous reference: G, F, K, L against CSR ----------------------
# Reported for completeness and as the straw man the report then dismantles.
# Border correction is used for K and L: it is cheap, unbiased, and at these
# sample sizes the efficiency gain from Ripley's isotropic correction does not
# justify its cost.
t0 <- Sys.time()
env_G <- envelope(X, Gest, r = r_GF, correction = "km",   nsim = NSIM,
                  savefuns = FALSE, verbose = FALSE)
env_F <- envelope(X, Fest, r = r_GF, correction = "km",   nsim = NSIM,
                  savefuns = FALSE, verbose = FALSE)
env_L <- envelope(X, Lest, r = r_K,  correction = "border", nsim = NSIM,
                  savefuns = FALSE, verbose = FALSE)
message(sprintf("CSR envelopes done in %.1f min", as.numeric(difftime(Sys.time(), t0, units="mins"))))

# --- 2. pair correlation: the sensor's own footprint ------------------------
# g(r) is a density, not a cumulative quantity, so unlike K it localises scale.
# It is the statistic that exposes the 375 m pixel floor: no two detections can
# be arbitrarily close, because the sensor cannot resolve them.
pcf_obs <- pcf(X, r = r_pcf, divisor = "d", correction = "translate")

# --- 3. the honest test: inhomogeneous second-order statistics --------------
# lambda-hat is re-estimated for each pattern, data and simulated alike, with a
# fixed bandwidth. Re-selecting the bandwidth inside the loop would be more
# principled still, but the cost is prohibitive and the fixed value is applied
# identically to data and simulations, so no asymmetry is introduced.
lam_data <- density(X, sigma = sigma, edge = TRUE, diggle = TRUE, positive = TRUE)

# `envelope` forwards its own `correction` argument to the summary function, so
# the setting is injected by overwriting the forwarded value rather than by
# passing it positionally, which would collide.
Linhom_reest <- function(Y, ...) {
  lam <- density(Y, sigma = sigma, edge = TRUE, diggle = TRUE, positive = TRUE)
  a <- list(...); a$correction <- "border"
  do.call(Linhom, c(list(Y, lambda = lam), a))
}
pcfinhom_reest <- function(Y, ...) {
  lam <- density(Y, sigma = sigma, edge = TRUE, diggle = TRUE, positive = TRUE)
  a <- list(...); a$correction <- "translate"; a$divisor <- "d"
  do.call(pcfinhom, c(list(Y, lambda = lam), a))
}

t0 <- Sys.time()
env_Linhom <- envelope(X, Linhom_reest, r = r_K, nsim = NSIM,
                       simulate = expression(rpoispp(lam_data)),
                       savefuns = FALSE, verbose = FALSE)
message(sprintf("Linhom envelope done in %.1f min", as.numeric(difftime(Sys.time(), t0, units="mins"))))
t0 <- Sys.time()
env_pcfinhom <- envelope(X, pcfinhom_reest, r = r_pcf, nsim = NSIM,
                         simulate = expression(rpoispp(lam_data)),
                         savefuns = FALSE, verbose = FALSE)
message(sprintf("pcfinhom envelope done in %.1f min", as.numeric(difftime(Sys.time(), t0, units="mins"))))

# --- 4. a global test, not a pointwise band --------------------------------
# Pointwise envelopes are read at every r, which inflates the type-I error.
# The maximum absolute deviation test gives one honest p-value for the whole
# curve, so the report can make a claim of significance that survives scrutiny.
# correction is pinned to "border" to match the envelopes above. Left unset,
# Lest computes border, Ripley and translation corrections on every one of the
# simulated patterns, which is both far slower and inconsistent with the
# correction actually reported.
Lest_border <- function(Y, ...) { a <- list(...); a$correction <- "border"
                                  do.call(Lest, c(list(Y), a)) }
mad_L      <- mad.test(X, Lest_border, r = r_K, nsim = 99, verbose = FALSE)
mad_Linhom <- mad.test(X, Linhom_reest, r = r_K, nsim = 99,
                       simulate = expression(rpoispp(lam_data)), verbose = FALSE)
dclf_Linhom <- dclf.test(X, Linhom_reest, r = r_K, nsim = 99,
                         simulate = expression(rpoispp(lam_data)), verbose = FALSE)

# --- 5. sensitivity to the event definition --------------------------------
# The same inhomogeneous statistic on the un-thinned pattern. If the conclusion
# flips, the finding is an artefact of the event definition and must be reported
# as such.
lam_all <- density(XA, sigma = sigma, edge = TRUE, diggle = TRUE, positive = TRUE)
Linhom_all <- Linhom(XA, lambda = lam_all, r = r_K, correction = "border")
Linhom_cell <- Linhom(X, lambda = lam_data, r = r_K, correction = "border")

# fv and envelope objects carry the closures that built them, and those closures
# capture the calling environment, which here includes other multi-megabyte caches.
# Storing plain data frames keeps every number and drops hundreds of megabytes.
slim_fv <- function(x) if (inherits(x, c("fv", "envelope"))) as.data.frame(x) else x
saveRDS(lapply(list(
  r_K = r_K, r_GF = r_GF, r_pcf = r_pcf, nsim = NSIM, sigma = sigma,
  RMAX_K = RMAX_K, RMAX_GF = RMAX_GF, RMAX_PCF = RMAX_PCF,
  env_G = env_G, env_F = env_F, env_L = env_L, pcf_obs = pcf_obs,
  lam_data = lam_data,
  env_Linhom = env_Linhom, env_pcfinhom = env_pcfinhom,
  mad_L = mad_L, mad_Linhom = mad_Linhom, dclf_Linhom = dclf_Linhom,
  Linhom_all = Linhom_all, Linhom_cell = Linhom_cell,
  n_cell = npoints(X), n_all = npoints(XA)
), slim_fv), file.path(DIR_DERIVED, "secondorder.rds"))
message("second-order analysis cached")
