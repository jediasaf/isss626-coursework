# ---------------------------------------------------------------------------
# 13 - Two corrections to R/12
#
# (a) The land-cover distance maps were undefined over 98% of the window. distmap
#     returns distances on the enclosing frame of its argument, and the argument
#     was the union of the landuse polygons, whose frame is their own bounding
#     box. Evaluating distfun over the regency window instead fixes it.
#
# (b) The non-parametric bootstrap in R/12 resampled points with replacement.
#     That is not a valid bootstrap for a point process: it destroys the spatial
#     dependence being estimated and creates coincident points, which distorts
#     the pair correlation function at exactly the short lags the fit relies on.
#     The symptom was clear, the point estimate fell outside the resulting
#     interval. This replaces it with a parametric bootstrap: simulate from the
#     fitted inhomogeneous Thomas process, refit each simulation by the same
#     minimum contrast, and take the spread of the refits.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(purrr); library(glue)
  library(spatstat.geom); library(spatstat.explore); library(spatstat.model)
  library(spatstat.random)
})
source("R/00_config.R")
set.seed(SEED)

P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
FO <- readRDS(file.path(DIR_DERIVED, "firstorder.rds"))
CM <- readRDS(file.path(DIR_DERIVED, "covariates_models.rds"))
TH <- readRDS(file.path(DIR_DERIVED, "thomas.rds"))
EX <- readRDS(file.path(DIR_DERIVED, "extended.rds"))
L12 <- readRDS(file.path(DIR_DERIVED, "landcover_ci.rds"))
X <- unmark(P$ppp_cell); win <- P$win
lam <- FO$kde$main

# --- (a) land-cover distances, over the whole window -----------------------
lu <- readRDS(file.path(DIR_DERIVED, "landcover.rds")) |>
  st_transform(STUDY_CRS) |> st_make_valid()
lu <- lu[lengths(st_intersects(lu, P$win_p)) > 0, ]
GROUPS <- list(forest = "forest",
               plantation = c("orchard", "farmland", "farmyard", "vineyard"))

dist_to_group <- function(classes) {
  g <- lu |> filter(fclass %in% classes)
  if (!nrow(g)) return(NULL)
  w <- as.owin(st_union(st_geometry(g)))
  # distfun evaluated on the analysis window, so the covariate is defined
  # everywhere a quadrature point can fall
  as.im(distfun(w), W = win, dimyx = c(600, 300))
}
D_cover <- compact(map(GROUPS, dist_to_group))
cover_na <- imap_dfr(D_cover, function(d, nm)
  tibble(group = nm, undefined_share = mean(is.na(as.vector(d$v)[!is.na(as.vector(as.im(win, dimyx=dim(d$v))$v))]))))
message("land-cover covariates rebuilt; checking coverage")
print(as.data.frame(cover_na))

rho_cover <- imap(D_cover, function(d, nm) as.data.frame(rhohat(X, d, confidence = 0.95)))

dr <- CM$D_road; dw <- CM$D_water
Droad_km <- eval.im(dr / 1000); Dwater_km <- eval.im(dw / 1000)
covs <- list(Droad_km = Droad_km, Dwater_km = Dwater_km)
for (nm in names(D_cover)) { di <- D_cover[[nm]]; covs[[paste0("D", nm, "_km")]] <- eval.im(di / 1000) }

fitm <- function(f) ppm(as.formula(f), data = covs)
mods <- list(
  null   = ppm(X ~ 1),
  road   = fitm("X ~ Droad_km"),
  access = fitm("X ~ Droad_km + Dwater_km"),
  cover  = fitm(paste("X ~", paste(grep("^D(forest|plantation)", names(covs), value = TRUE), collapse = " + "))),
  all    = fitm(paste("X ~", paste(names(covs), collapse = " + ")))
)
aic2 <- tibble(model = names(mods), AIC = sapply(mods, AIC),
               npar = sapply(mods, function(m) length(coef(m)))) |>
  mutate(dAIC = AIC - min(AIC)) |> arrange(AIC)
print(as.data.frame(aic2))
best2 <- tryCatch(coef(summary(mods[[aic2$model[1]]])), error = function(e) NULL)

# --- (b) parametric bootstrap on the cluster scale -------------------------
RMIN <- TH$RMIN; RMAX <- TH$RMAX; Q <- TH$Q
tp <- TH$thomas_tab |> filter(key == "trend_removed")
A_hat <- tp$A; sigma_hat <- tp$sigma_m
kappa_hat <- 1 / (4 * pi * A_hat * sigma_hat^2)     # per square metre
# Offspring intensity image so that kappa * mu(u) reproduces lambda-hat(u).
mu_im <- eval.im(lam / kappa_hat)

thomas_g <- function(r, A, sigma) 1 + A * exp(-r^2 / (4 * sigma^2))
fit_once <- function(r, ghat) {
  ok <- is.finite(r) & is.finite(ghat) & r >= RMIN & r <= RMAX
  r <- r[ok]; ghat <- ghat[ok]
  if (length(r) < 20) return(NA_real_)
  obj <- function(par) sum((ghat^Q - thomas_g(r, exp(par[1]), exp(par[2]))^Q)^2)
  best <- NULL
  for (st in list(c(log(2), log(700)), c(log(5), log(1500)), c(log(1), log(300)))) {
    o <- try(optim(st, obj, method = "Nelder-Mead",
                   control = list(maxit = 2000, reltol = 1e-10)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) NA_real_ else exp(best$par[2])
}

NB <- as.integer(Sys.getenv("NPARBOOT", 199L))
rgrid <- EX$r_ext[EX$r_ext <= RMAX]
message(glue("parametric bootstrap from the fitted Thomas process, {NB} simulations"))
sig_par <- numeric(NB)
for (b in seq_len(NB)) {
  Yb <- try(rThomas(kappa = kappa_hat, scale = sigma_hat, mu = mu_im, win = win),
            silent = TRUE)
  if (inherits(Yb, "try-error") || npoints(Yb) < 500) { sig_par[b] <- NA; next }
  lamb <- density(Yb, sigma = FO$sigma_main, edge = TRUE, diggle = TRUE, positive = TRUE)
  gb <- try(pcfinhom(Yb, lambda = lamb, r = rgrid, divisor = "d",
                     correction = "translate"), silent = TRUE)
  if (inherits(gb, "try-error")) { sig_par[b] <- NA; next }
  gb <- as.data.frame(gb)
  sig_par[b] <- fit_once(gb$r, gb$trans)
  if (b %% 50 == 0) message(glue("  {b}/{NB}"))
}
sig_par <- sig_par[is.finite(sig_par)]
q <- quantile(sig_par, c(0.025, 0.5, 0.975))
par_boot <- tibble(
  key = "par",
  method = "parametric, from the fitted Thomas process",
  point_sigma_m = sigma_hat, median_m = unname(q[2]),
  lo95_m = unname(q[1]), hi95_m = unname(q[3]),
  r95_point_m = 2.4477 * sigma_hat,
  r95_lo_m = 2.4477 * unname(q[1]), r95_hi_m = 2.4477 * unname(q[3]),
  n_ok = length(sig_par), n_sim = NB,
  mean_n_sim = NA_real_)
print(as.data.frame(par_boot |> mutate(across(where(is.numeric), ~ signif(.x, 4)))))

# the invalid non-parametric interval, kept for the comparison in the report
np_boot <- L12$boot |> transmute(
  key = "np",
  method = "non-parametric, resampling points with replacement",
  point_sigma_m, median_m = boot_median_m, lo95_m, hi95_m,
  r95_point_m = 2.4477 * point_sigma_m, r95_lo_m, r95_hi_m,
  n_ok, n_sim = n_boot, mean_n_sim = NA_real_)

saveRDS(list(aic2 = aic2, best_coef2 = best2, rho_cover = rho_cover,
             cover_share = L12$cover_share, cover_tab = L12$cover_tab,
             boot_compare = bind_rows(par_boot, np_boot),
             sigma_par = sig_par, sigma_np = L12$sigma_boot,
             kappa_hat = kappa_hat),
        file.path(DIR_DERIVED, "cover_boot.rds"))
message("corrected land cover and parametric bootstrap cached")
