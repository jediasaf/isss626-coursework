# ---------------------------------------------------------------------------
# 06 - Value-added: covariates, fitted models, and a cluster-scale estimate
#
# The core tasks establish that fire is clustered, where, and over what distance.
# None of that says what the clustering is organised around or how large a
# cluster is, which is the form a planner can use. Three extensions:
#
#  A  Covariate dependence. Peat fire in the former Mega Rice Project area is
#     overwhelmingly anthropogenic in ignition and made possible by drainage.
#     Distance to the road network proxies access and ignition opportunity;
#     distance to mapped canals and waterways proxies drainage-driven peat
#     desiccation. rhohat estimates intensity as a non-parametric function of
#     each, imposing no functional form.
#
#  B  Fitted point-process models, compared by AIC, with the residual K function
#     used to show what a Poisson model still fails to capture.
#
#  C  Cluster scale, from a Thomas process fitted by minimum contrast. Unlike a
#     significance test this returns interpretable quantities: parent intensity,
#     mean events per cluster, and a cluster radius.
#
# Covariate geometry comes from a Geofabrik OpenStreetMap extract rather than a
# live Overpass query. Overpass repeatedly refused a bounding box this large, and
# a pinned dated extract is in any case the reproducible choice: a live query
# returns different data every time it is run.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(spatstat.geom); library(spatstat.explore)
  library(spatstat.model); library(spatstat.random)
})
source("R/00_config.R")
set.seed(SEED)
P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
FO <- readRDS(file.path(DIR_DERIVED, "firstorder.rds"))
# ppm, rhohat, kppm and scan.test operate on unmarked or multitype patterns only.
# The analysis marks (FRP, confidence, platform, date) are carried for other
# purposes and are stripped here rather than being silently ignored with a warning.
X  <- unmark(P$ppp_cell)
win <- P$win

# --- A. OpenStreetMap covariates from the Geofabrik Kalimantan extract ------
zipf <- file.path(DIR_RAW, "kalimantan-latest-free.gpkg.zip")
osm_file <- file.path(DIR_DERIVED, "osm_covariates.rds")

if (!file.exists(osm_file)) {
  stopifnot(file.exists(zipf))
  gpkg <- grep("\\.gpkg$", unzip(zipf, list = TRUE)$Name, value = TRUE)[1]
  vsi  <- paste0("/vsizip/", normalizePath(zipf), "/", gpkg)
  lyrs <- st_layers(vsi)$name
  message("layers: ", paste(lyrs, collapse = ", "))

  # A generous buffer around the regency, so distance-to-nearest is correct for
  # points near the boundary: the nearest road to a southern cell may lie outside
  # the regency entirely, and clipping first would bias every such distance upward.
  buf <- st_buffer(st_transform(P$win_sf, STUDY_CRS), 15000) |> st_transform(4326)
  wkt <- st_as_text(st_as_sfc(st_bbox(buf)))

  pick <- function(pattern) {
    l <- grep(pattern, lyrs, value = TRUE)[1]
    if (is.na(l)) return(NULL)
    message("reading ", l)
    st_read(vsi, layer = l, wkt_filter = wkt, quiet = TRUE)
  }
  roads <- pick("roads")
  water <- pick("waterways")
  saveRDS(list(roads = roads, water = water), osm_file)
}
osm <- readRDS(osm_file)
message(sprintf("roads features: %s | waterway features: %s",
                format(NROW(osm$roads), big.mark = ","),
                format(NROW(osm$water), big.mark = ",")))
road_types  <- if (!is.null(osm$roads)) sort(table(osm$roads$fclass), decreasing = TRUE) else NULL
water_types <- if (!is.null(osm$water)) sort(table(osm$water$fclass), decreasing = TRUE) else NULL

# The distance map must be computed on a frame larger than the analysis window,
# then cropped, for the boundary reason above.
frame_p <- st_buffer(st_transform(P$win_sf, STUDY_CRS), 15000)
frame_w <- as.owin(st_geometry(frame_p))

to_psp <- function(lines, w) {
  if (is.null(lines) || NROW(lines) == 0) return(NULL)
  g <- st_geometry(st_transform(lines, STUDY_CRS))
  g <- suppressWarnings(st_intersection(g, st_geometry(frame_p)))
  g <- g[!st_is_empty(g)]
  g <- suppressWarnings(st_cast(g, "LINESTRING", warn = FALSE))
  if (length(g) == 0) return(NULL)
  segs <- do.call(rbind, lapply(seq_along(g), function(k) {
    m <- st_coordinates(g[[k]])
    if (is.null(dim(m)) || nrow(m) < 2) return(NULL)
    m <- m[, 1:2, drop = FALSE]
    cbind(m[-nrow(m), , drop = FALSE], m[-1, , drop = FALSE])
  }))
  if (is.null(segs)) return(NULL)
  psp(segs[,1], segs[,2], segs[,3], segs[,4], window = w, check = FALSE)
}
road_psp  <- to_psp(osm$roads, frame_w)
water_psp <- to_psp(osm$water, frame_w)
message(sprintf("road segments: %s | waterway segments: %s",
                ifelse(is.null(road_psp),  "0", format(road_psp$n,  big.mark = ",")),
                ifelse(is.null(water_psp), "0", format(water_psp$n, big.mark = ","))))

mk_dist <- function(p) {
  if (is.null(p)) return(NULL)
  d <- distmap(p, dimyx = c(600, 300))
  d[win, drop = FALSE]                      # crop to the analysis window
}
D_road  <- mk_dist(road_psp)
D_water <- mk_dist(water_psp)

rho_road  <- if (!is.null(D_road))  rhohat(X, D_road,  confidence = 0.95) else NULL
rho_water <- if (!is.null(D_water)) rhohat(X, D_water, confidence = 0.95) else NULL

# --- B. fitted Poisson models ----------------------------------------------
mods <- list(); mods$null <- ppm(X ~ 1)
if (!is.null(D_road) && !is.null(D_water)) {
  Droad_km  <- eval.im(D_road  / 1000)
  Dwater_km <- eval.im(D_water / 1000)
  mods$road  <- ppm(X ~ Droad_km)
  mods$water <- ppm(X ~ Dwater_km)
  mods$both  <- ppm(X ~ Droad_km + Dwater_km)
  # A log term lets a steep near-field effect flatten, which a linear term in
  # distance cannot represent and which rhohat is expected to show.
  mods$log   <- ppm(X ~ log(Droad_km + 0.1) + log(Dwater_km + 0.1))
}
aic <- sapply(mods, AIC)
aic_tab <- tibble(model = names(aic), AIC = as.numeric(aic),
                  dAIC = as.numeric(aic) - min(as.numeric(aic)),
                  npar = sapply(mods, function(m) length(coef(m)))) |> arrange(AIC)
print(as.data.frame(aic_tab))
best_cov <- mods[[aic_tab$model[1]]]

res_K <- tryCatch(Kres(best_cov, correction = "border"), error = function(e) NULL)

# --- C. Thomas cluster process ---------------------------------------------
fit_thomas <- tryCatch(
  kppm(X ~ 1, clusters = "Thomas", statistic = "pcf",
       statargs = list(divisor = "d"), rmax = 3000),
  error = function(e) { message("kppm Thomas failed: ", conditionMessage(e)); NULL })
fit_thomas_trend <- if (!is.null(D_road) && !is.null(D_water)) tryCatch(
  kppm(X ~ Droad_km + Dwater_km, clusters = "Thomas", statistic = "pcf",
       statargs = list(divisor = "d"), rmax = 3000),
  error = function(e) NULL) else NULL

thomas_par <- if (!is.null(fit_thomas)) {
  p <- fit_thomas$clustpar
  tibble(kappa_per_km2 = unname(p["kappa"]) * 1e6,
         sigma_m = unname(p["scale"]),
         mu_per_cluster = unname(fit_thomas$mu),
         # For a bivariate Gaussian offspring kernel, ~95% of offspring lie
         # within 2.4477 sigma of the parent. A usable operational radius.
         r95_m = 2.4477 * unname(p["scale"]))
} else NULL
if (!is.null(thomas_par)) print(as.data.frame(thomas_par))

# --- D. spatial scan statistic ---------------------------------------------
scan_res <- tryCatch(
  scan.test(X, r = c(2000, 5000, 10000), method = "poisson", nsim = 99, verbose = FALSE),
  error = function(e) { message("scan.test failed: ", conditionMessage(e)); NULL })
if (!is.null(scan_res)) print(scan_res)

# Only what the report and slides read is cached. Retaining the fitted ppm objects
# and the 151,558-segment line patterns produced a 387 MB file, which would slow
# every render and is pointless: the coefficients, AIC table and residual summary
# are the outputs, and the models are re-fittable by re-running this script.
thomas_summary <- if (!is.null(fit_thomas)) capture.output(print(fit_thomas)) else NULL
saveRDS(list(
  road_types = road_types, water_types = water_types,
  n_roads = NROW(osm$roads), n_water = NROW(osm$water),
  n_road_seg = if (is.null(road_psp)) 0L else road_psp$n,
  n_water_seg = if (is.null(water_psp)) 0L else water_psp$n,
  D_road = D_road, D_water = D_water,
  rho_road = rho_road, rho_water = rho_water,
  aic_tab = aic_tab, best_model = aic_tab$model[1],
  best_coef = tryCatch(coef(summary(best_cov)), error = function(e) NULL),
  res_K = res_K, thomas_par = thomas_par, thomas_summary = thomas_summary,
  scan_res = scan_res
), file.path(DIR_DERIVED, "covariates_models.rds"))
message("covariate and model analysis cached")
