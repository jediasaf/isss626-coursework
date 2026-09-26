# ---------------------------------------------------------------------------
# 12 - Land cover as a third covariate, and a confidence interval on the
#      cluster scale
#
# Two gaps left over from the core analysis.
#
# The covariate set was thin: distance to road and distance to mapped waterway,
# no land cover at all, in a report whose entire framing is peatland. The
# Geofabrik extract carries OSM landuse polygons, which are coarse but real, and
# they are already on disk.
#
# The cluster scale was quoted as a point estimate. Minimum contrast gives no
# standard error, so the number was unqualified. A bootstrap over the pattern
# supplies one.
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
X  <- unmark(P$ppp_cell)
win <- P$win

# --- A. land cover ----------------------------------------------------------
lc_file <- file.path(DIR_DERIVED, "landcover.rds")
if (!file.exists(lc_file)) {
  z <- file.path(DIR_RAW, "kalimantan-latest-free.gpkg.zip")
  g <- grep("[.]gpkg$", unzip(z, list = TRUE)$Name, value = TRUE)[1]
  v <- paste0("/vsizip/", normalizePath(z), "/", g)
  wkt <- st_as_text(st_as_sfc(st_bbox(st_transform(P$win_sf, 4326))))
  lu <- st_read(v, layer = "gis_osm_landuse_a_free", wkt_filter = wkt, quiet = TRUE)
  saveRDS(lu, lc_file)
}
lu <- readRDS(lc_file)
message(glue("OSM landuse polygons in the regency bbox: {nrow(lu)}"))

lu_p <- lu |> st_transform(STUDY_CRS) |> st_make_valid()
lu_p <- lu_p[lengths(st_intersects(lu_p, P$win_p)) > 0, ]
cover_tab <- lu_p |> st_drop_geometry() |> count(fclass, sort = TRUE)
print(as.data.frame(head(cover_tab, 12)))

# The classes that matter for this landscape. Everything else is lumped, because
# OSM coverage of the rest is too patchy to model separately.
GROUPS <- list(forest = c("forest"),
               plantation = c("orchard", "farmland", "farmyard", "vineyard"),
               wetland = c("wetland", "marsh", "scrub", "meadow", "grass"))

as_mask <- function(classes) {
  g <- lu_p |> filter(fclass %in% classes)
  if (!nrow(g)) return(NULL)
  u <- st_union(st_geometry(g))
  as.owin(st_intersection(u, st_geometry(P$win_p)))
}
masks <- compact(map(GROUPS, as_mask))
message(glue("land-cover groups usable: {paste(names(masks), collapse=', ')}"))

# Distance to the nearest polygon of each group, as a covariate image. Distance
# rather than a categorical indicator, because OSM polygons do not tile the
# regency and a categorical covariate would have a huge undefined class.
D_cover <- imap(masks, function(w, nm) {
  d <- distmap(w, dimyx = c(600, 300))
  d[win, drop = FALSE]
})
rho_cover <- imap(D_cover, function(d, nm) rhohat(X, d, confidence = 0.95))

# Fraction of the window each group covers, which sets how much any of this can
# possibly explain.
cover_share <- imap_dfr(masks, function(w, nm)
  tibble(group = nm, share_of_window = area.owin(w) / area.owin(win)))
print(as.data.frame(cover_share))

# --- model comparison, now with land cover ---------------------------------
# eval.im resolves images by symbol, so each image is bound to a plain name first.
dr <- CM$D_road; dw <- CM$D_water
Droad_km  <- eval.im(dr / 1000)
Dwater_km <- eval.im(dw / 1000)
covs <- list(Droad_km = Droad_km, Dwater_km = Dwater_km)
for (nm in names(D_cover)) {
  di <- D_cover[[nm]]
  covs[[paste0("D", nm, "_km")]] <- eval.im(di / 1000)
}

fit <- function(f) ppm(as.formula(f), data = covs)
mods <- list(
  null  = ppm(X ~ 1),
  road  = fit("X ~ Droad_km"),
  both  = fit("X ~ Droad_km + Dwater_km"),
  cover = fit(paste("X ~", paste(setdiff(names(covs), c("Droad_km","Dwater_km")), collapse = " + "))),
  all   = fit(paste("X ~", paste(names(covs), collapse = " + ")))
)
aic2 <- tibble(model = names(mods), AIC = sapply(mods, AIC),
               npar = sapply(mods, function(m) length(coef(m)))) |>
  mutate(dAIC = AIC - min(AIC)) |> arrange(AIC)
print(as.data.frame(aic2))
best_coef2 <- tryCatch(coef(summary(mods[[aic2$model[1]]])), error = function(e) NULL)

# --- B. bootstrap interval on the cluster scale ----------------------------
# The Thomas fit in R/07 minimises a contrast between the observed inhomogeneous
# pcf and the model pcf. There is no analytic standard error for that, so the
# interval comes from refitting on bootstrap resamples of the pattern. Points are
# resampled with replacement inside the same window; lambda-hat is recomputed each
# time so the trend removal is resampled along with the interaction.
EX <- readRDS(file.path(DIR_DERIVED, "extended.rds"))
RMIN <- TH$RMIN; RMAX <- TH$RMAX; Q <- TH$Q
thomas_g <- function(r, A, sigma) 1 + A * exp(-r^2 / (4 * sigma^2))

fit_mc_once <- function(r, ghat) {
  ok <- is.finite(r) & is.finite(ghat) & r >= RMIN & r <= RMAX
  r <- r[ok]; ghat <- ghat[ok]
  if (length(r) < 20) return(c(NA, NA))
  obj <- function(par) sum((ghat^Q - thomas_g(r, exp(par[1]), exp(par[2]))^Q)^2)
  best <- NULL
  for (st in list(c(log(2), log(700)), c(log(5), log(1500)), c(log(1), log(300)))) {
    o <- try(optim(st, obj, method = "Nelder-Mead",
                   control = list(maxit = 2000, reltol = 1e-10)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(c(NA, NA))
  c(exp(best$par[1]), exp(best$par[2]))
}

NBOOT <- as.integer(Sys.getenv("NBOOT_THOMAS", 199L))
sigma_boot <- numeric(NBOOT)
n <- npoints(X); co <- coords(X)
message(glue("bootstrapping the Thomas fit, {NBOOT} resamples"))
for (b in seq_len(NBOOT)) {
  idx <- sample.int(n, n, replace = TRUE)
  Yb <- ppp(co$x[idx], co$y[idx], window = win, checkdup = FALSE)
  lamb <- density(Yb, sigma = FO$sigma_main, edge = TRUE, diggle = TRUE, positive = TRUE)
  gb <- try(pcfinhom(Yb, lambda = lamb, r = EX$r_ext[EX$r_ext <= RMAX],
                     divisor = "d", correction = "translate"), silent = TRUE)
  if (inherits(gb, "try-error")) { sigma_boot[b] <- NA; next }
  gb <- as.data.frame(gb)
  sigma_boot[b] <- fit_mc_once(gb$r, gb$trans)[2]
  if (b %% 50 == 0) message(glue("  {b}/{NBOOT}"))
}
sigma_boot <- sigma_boot[is.finite(sigma_boot)]
point_est <- TH$thomas_tab$sigma_m[TH$thomas_tab$key == "trend_removed"]
ci <- quantile(sigma_boot, c(0.025, 0.5, 0.975))
boot_summary <- tibble(
  point_sigma_m = point_est,
  boot_median_m = unname(ci[2]),
  lo95_m = unname(ci[1]), hi95_m = unname(ci[3]),
  r95_lo_m = 2.4477 * unname(ci[1]), r95_hi_m = 2.4477 * unname(ci[3]),
  n_ok = length(sigma_boot), n_boot = NBOOT)
print(as.data.frame(boot_summary |> mutate(across(everything(), ~ signif(.x, 4)))))

saveRDS(list(cover_tab = cover_tab, cover_share = cover_share,
             rho_cover = map(rho_cover, as.data.frame),
             aic2 = aic2, best_coef2 = best_coef2,
             sigma_boot = sigma_boot, boot = boot_summary,
             groups = names(masks)),
        file.path(DIR_DERIVED, "landcover_ci.rds"))
message("land cover and bootstrap cached")
