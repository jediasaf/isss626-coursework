# ---------------------------------------------------------------------------
# 03 - First-order properties: how intensity varies in space and in time
#
# Method choices and why:
#  * Quadrat counts give a blunt, interpretable test of a homogeneous Poisson
#    process. They are reported because they are transparent, and criticised
#    because the verdict depends on an arbitrary quadrat size; the sensitivity
#    across sizes is shown rather than a single favourable choice.
#  * Kernel density estimation is the substantive first-order tool. The
#    bandwidth is the decision that matters, so four selectors are compared and
#    the choice is argued, not defaulted.
#  * Edge correction is applied: the regency is long and narrow (86 x 213 km),
#    so an uncorrected estimate would understate intensity along the boundary.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(spatstat.geom); library(spatstat.explore)
})
source("R/00_config.R")
set.seed(SEED)
P <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
X <- P$ppp_cell          # primary spatial pattern: one point per 375 m cell
XA <- P$ppp_all          # sensitivity: every retained detection
message(sprintf("ppp_cell n=%d  ppp_all n=%d  area=%.0f km2",
                npoints(X), npoints(XA), P$area_km2))

# --- 1. quadrat counts, across quadrat sizes -------------------------------
# The window is roughly 1:2.5, so quadrats are set to be near-square in ground
# units rather than an equal nx/ny split, which would produce elongated cells.
quad_grid <- tibble::tibble(nx = c(3L, 5L, 8L, 12L, 16L)) |>
  mutate(ny = pmax(1L, round(nx * 213.3 / 86.1)))
quad_tests <- lapply(seq_len(nrow(quad_grid)), function(i) {
  qt <- quadrat.test(X, nx = quad_grid$nx[i], ny = quad_grid$ny[i])
  list(nx = quad_grid$nx[i], ny = quad_grid$ny[i], test = qt,
       stat = unname(qt$statistic), df = unname(qt$parameter),
       p = qt$p.value, n_quadrats = qt$parameter + 1)
})
quad_summary <- do.call(rbind, lapply(quad_tests, function(z)
  data.frame(nx = z$nx, ny = z$ny, n_quadrats = z$n_quadrats,
             X2 = z$stat, df = z$df, p_value = z$p)))
print(quad_summary)

# Variance-to-mean ratio, a scale-explicit description that does not pretend to
# be a hypothesis test. VMR = 1 under CSR, > 1 indicates over-dispersion.
vmr <- sapply(seq_len(nrow(quad_grid)), function(i) {
  cnt <- as.vector(quadratcount(X, nx = quad_grid$nx[i], ny = quad_grid$ny[i]))
  var(cnt) / mean(cnt)
})
quad_summary$vmr <- round(vmr, 1)

# --- 2. bandwidth selection ------------------------------------------------
# bw.diggle    minimises a MSE criterion for the intensity of a Cox process;
#              tends to select a small bandwidth, good for detecting fine
#              structure, prone to noise.
# bw.ppl       likelihood cross-validation; usually larger, better when the
#              target is the broad intensity trend.
# bw.CvL       Cronie and van Lieshout's criterion, tuned so the estimated
#              intensity integrates to the observed count.
# bw.scott     a fast rule of thumb, no optimisation; included as a reference.
bw_tab <- tibble::tibble(
  selector = c("bw.diggle", "bw.ppl", "bw.CvL", "bw.scott (x)", "bw.scott (y)"),
  sigma_m = c(as.numeric(bw.diggle(X)), as.numeric(bw.ppl(X)),
              as.numeric(bw.CvL(X)), as.numeric(bw.scott(X))[1],
              as.numeric(bw.scott(X))[2])
) |> mutate(sigma_km = round(sigma_m / 1000, 2))
print(bw_tab)

# --- the bandwidth actually adopted, and why not a cross-validated one -------
# Both cross-validation selectors return sub-kilometre bandwidths, and
# density.ppp emits a numerical-underflow warning at that scale over a window
# this large. That is not a defect in the selectors; it is the documented
# behaviour of likelihood and MSE cross-validation when a pattern is strongly
# clustered at small scales. The criteria reward a surface that reproduces the
# local bursts, so they interpret clustering as fine-grained intensity structure.
#
# Adopting such a bandwidth here would be self-defeating. Section 7 tests whether
# residual interaction exists once first-order intensity is accounted for, using
# this same estimate as the reference intensity. If lambda-hat already absorbs
# the sub-kilometre clustering, the second-order test is asking whether the
# clustering exists after removing it, and the answer is rigged.
#
# The adopted bandwidth is therefore fixed at 4 km, on a criterion stated in
# advance: it must exceed the interaction range so that trend and interaction are
# estimated by separate parts of the analysis rather than competing for the same
# signal. 4 km sits above the ~2 km interaction range found in Section 7, inside
# the interval spanned by bw.scott's two components, and below bw.CvL. All
# selectors are reported, and the comparison panel shows what each would produce.
sigma_main   <- 4000
sigma_ppl    <- as.numeric(bw.ppl(X))
sigma_diggle <- as.numeric(bw.diggle(X))
sigma_CvL    <- as.numeric(bw.CvL(X))
sigma_scott  <- as.numeric(bw.scott(X))
spatstat.options(npixel = c(256, 512))   # match the window's 1:2.5 aspect

kde <- list(
  ppl     = suppressWarnings(density(X, sigma = sigma_ppl, edge = TRUE,
                                     diggle = TRUE, positive = TRUE)),
  main    = density(X, sigma = sigma_main, edge = TRUE, diggle = TRUE, positive = TRUE),
  scott   = density(X, sigma = sigma_scott[1], edge = TRUE, diggle = TRUE, positive = TRUE),
  CvL     = density(X, sigma = sigma_CvL, edge = TRUE, diggle = TRUE, positive = TRUE)
)
# Adaptive estimator: bandwidth varies with local point density, so sparse
# areas are not over-smoothed and dense areas are not blurred. Reported as a
# check on whether the fixed-bandwidth picture is an artefact of one scale.
kde$adaptive <- adaptive.density(X, method = "voronoi", f = 0.1, nrep = 20)

# im$v is per square metre; 1e9 converts to per 1,000 km2, the reported unit.
kde_range <- sapply(kde, function(im) range(as.vector(im$v), na.rm = TRUE) * 1e9)
message("KDE intensity range, burned cells per 1,000 km2:"); print(round(kde_range))

# --- 3. intensity weighted by radiative power ------------------------------
# A count surface treats a smouldering peat pixel and an intense flaming front
# as equal. FRP (MW) is the sensor's estimate of radiative output, so a
# FRP-weighted surface answers a different and more operationally relevant
# question: where is energy being released, not merely where are detections.
Xfrp <- X
marks(Xfrp) <- as.numeric(marks(X)$frp)
Xfrp <- Xfrp[is.finite(marks(Xfrp)) & marks(Xfrp) > 0]
frp_smooth <- Smooth(Xfrp, sigma = sigma_main, edge = TRUE)
frp_total  <- density(Xfrp, sigma = sigma_main, weights = marks(Xfrp),
                      edge = TRUE, diggle = TRUE, positive = TRUE)

# --- 4. temporal first-order structure ------------------------------------
daily <- P$daily_daycell |>
  mutate(roll7 = as.numeric(stats::filter(n_events, rep(1/7, 7), sides = 2)),
         week  = as.integer(floor(as.integer(date_local - WINDOW_START) / 7)) + 1L,
         month = format(date_local, "%b"))
dn <- P$fires |> st_drop_geometry() |> count(date_local, daynight)
# Local solar time of the two overpasses, to show that "day" and "night" here
# are two fixed sampling windows, not a continuous diurnal record.
overpass <- P$fires |> st_drop_geometry() |>
  mutate(hr = as.numeric(format(ts_local, "%H")) +
              as.numeric(format(ts_local, "%M")) / 60) |>
  count(platform, daynight, hr_bin = round(hr * 2) / 2)

saveRDS(list(
  quad_summary = quad_summary, quad_tests = quad_tests, bw_tab = bw_tab,
  sigma_main = sigma_main, sigma_ppl = sigma_ppl, sigma_diggle = sigma_diggle,
  sigma_CvL = sigma_CvL, sigma_scott = sigma_scott,
  kde = kde, kde_range = kde_range,
  frp_smooth = frp_smooth, frp_total = frp_total,
  daily = daily, dn = dn, overpass = overpass
), file.path(DIR_DERIVED, "firstorder.rds"))
message("first-order analysis cached")
