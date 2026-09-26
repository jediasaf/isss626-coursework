# ---------------------------------------------------------------------------
# 07 - Cluster scale, by explicit minimum contrast
#
# A homogeneous Thomas process fitted to this pattern by kppm returns a cluster
# standard deviation of about 3.3 km, roughly 440 events per cluster, and only
# about 23 parent locations in the whole regency. That is not a cluster scale. It
# is the southern peat belt being re-described as a handful of enormous clusters,
# because a model with no trend term can only represent large-scale intensity
# variation by inventing very large clusters. This is the standard
# trend-versus-interaction identifiability problem in cluster processes.
#
# The fix is to remove the trend before estimating the cluster parameters. The
# inhomogeneous pair correlation function already does exactly that: it divides
# out lambda-hat(u). So the Thomas pcf is fitted directly to g_inhom.
#
# The fit is written out rather than delegated, because the objective function is
# the substance. For a Thomas process,
#
#     g(r) = 1 + (1 / (4 pi kappa sigma^2)) * exp(-r^2 / (4 sigma^2))
#
# and the minimum-contrast estimate minimises
#
#     sum over r in [rmin, rmax] of ( ghat(r)^q - g_theta(r)^q )^2
#
# with q = 1/4, the transformation spatstat uses by default, which stabilises the
# variance of the pcf estimate across r. kappa and the mean cluster size follow
# from the fitted amplitude and the observed mean intensity.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(spatstat.geom); library(spatstat.explore)
})
source("R/00_config.R")
set.seed(SEED)
P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
SO <- readRDS(file.path(DIR_DERIVED, "secondorder.rds"))
X  <- P$ppp_cell
area_m2   <- area.owin(P$win)
lambda_bar <- npoints(X) / area_m2

RMIN <- 200      # below this, g(r) is governed by the 375 m sensor floor, not by
                 # the cluster process, so those lags would bias the fit
RMAX <- 5000
Q    <- 0.25

thomas_g <- function(r, A, sigma) 1 + A * exp(-r^2 / (4 * sigma^2))

fit_mc <- function(r, ghat, label, key) {
  ok <- is.finite(r) & is.finite(ghat) & r >= RMIN & r <= RMAX
  r <- r[ok]; ghat <- ghat[ok]
  if (!length(r)) return(NULL)
  obj <- function(par) {
    A <- exp(par[1]); sigma <- exp(par[2])
    sum((ghat^Q - thomas_g(r, A, sigma)^Q)^2)
  }
  # Several starts, because minimum contrast on a pcf is not reliably unimodal.
  starts <- expand.grid(logA = log(c(0.5, 2, 5, 15)),
                        logs = log(c(300, 700, 1500, 3000)))
  best <- NULL
  for (k in seq_len(nrow(starts))) {
    o <- try(optim(c(starts$logA[k], starts$logs[k]), obj,
                   method = "Nelder-Mead",
                   control = list(maxit = 4000, reltol = 1e-12)), silent = TRUE)
    if (inherits(o, "try-error")) next
    if (is.null(best) || o$value < best$value) best <- o
  }
  if (is.null(best)) return(NULL)
  A <- exp(best$par[1]); sigma <- exp(best$par[2])
  kappa <- 1 / (4 * pi * A * sigma^2)
  tibble(key = key, spec = label, A = A, sigma_m = sigma,
         r95_m = 2.4477 * sigma,
         kappa_per_1000km2 = kappa * 1e9,
         n_parents_implied = kappa * area_m2,
         mu_per_cluster = lambda_bar / kappa,
         ss = best$value, n_r = length(r))
}

g_hom  <- as.data.frame(SO$pcf_obs)
g_inh  <- as.data.frame(SO$env_pcfinhom)

thomas_tab <- bind_rows(
  fit_mc(g_hom$r, g_hom$trans, "homogeneous pcf (trend NOT removed)", "no_trend"),
  fit_mc(g_inh$r, g_inh$obs,   "inhomogeneous pcf (trend removed)",   "trend_removed")
)
# The kppm result from R/06, retained for the contrast.
CMf <- file.path(DIR_DERIVED, "covariates_models.rds")
kppm_hom <- if (file.exists(CMf)) readRDS(CMf)$thomas_par else NULL

print(as.data.frame(thomas_tab |> mutate(across(where(is.numeric), ~ signif(.x, 4)))))
if (!is.null(kppm_hom)) {
  cat("\nkppm homogeneous fit from R/06 for comparison:\n")
  print(as.data.frame(kppm_hom |> mutate(across(everything(), ~ signif(.x, 4)))))
}

# Fitted curves, for plotting against the observed pcf.
r_plot <- seq(0, RMAX, length.out = 400)
curves <- bind_rows(lapply(seq_len(nrow(thomas_tab)), function(k)
  tibble(key = thomas_tab$key[k], spec = thomas_tab$spec[k], r = r_plot,
         g = thomas_g(r_plot, thomas_tab$A[k], thomas_tab$sigma_m[k]))))

# fv and envelope objects carry the closures that built them, and those closures
# capture the calling environment, which here includes other multi-megabyte caches.
# Storing plain data frames keeps every number and drops hundreds of megabytes.
slim_fv <- function(x) if (inherits(x, c("fv", "envelope"))) as.data.frame(x) else x
saveRDS(lapply(list(thomas_tab = thomas_tab, curves = curves, kppm_hom = kppm_hom,
             RMIN = RMIN, RMAX = RMAX, Q = Q, lambda_bar = lambda_bar),
        slim_fv), file.path(DIR_DERIVED, "thomas.rds"))
message("cluster-scale estimates cached")
