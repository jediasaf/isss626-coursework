# ---------------------------------------------------------------------------
# 11 - Testing assumption A7, isotropy
#
# I listed isotropy as an assumption and said anisotropy was likely, because peat
# fires spread with the wind and along canal lines. Then I never tested it. This
# script closes that gap.
#
# Why not the packaged tools. Both spatstat routines for this, Ksector and
# pairorient, are homogeneous-only: neither accepts a lambda argument. Applied
# directly here they would reject isotropy for two reasons that have nothing to
# do with fire:
#
#   * the window is 86 km wide and 213 km tall, so even a perfectly isotropic
#     process yields more north-south pairs than east-west ones simply because
#     there is more room in that direction
#   * intensity is concentrated in a south-southeast band, so pairs drawn from
#     that band inherit its elongation
#
# The test below therefore builds its own null by simulating from an
# inhomogeneous Poisson process with the same lambda-hat used in Section 7. That
# absorbs the window geometry and the first-order trend together, so what remains
# is directional structure in the interaction itself. It is the same move as
# replacing CSR with an inhomogeneous null, applied to direction instead of
# distance.
#
# Statistic. Orientations are axial, not directed: a pair at 10 degrees and a pair
# at 190 degrees describe the same alignment. Angles are therefore doubled before
# averaging, which is the standard device for axial data. The concentration
# measure is the mean resultant length of the doubled angles,
#
#     R = | mean( exp(2i*theta) ) |
#
# which is 0 for a perfectly isotropic set of orientations and 1 for perfect
# alignment. The preferred axis is 0.5 * arg( mean( exp(2i*theta) ) ).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(purrr); library(glue)
  library(spatstat.geom); library(spatstat.explore); library(spatstat.random)
})
source("R/00_config.R")
set.seed(SEED)

P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
FO <- readRDS(file.path(DIR_DERIVED, "firstorder.rds"))
X   <- unmark(P$ppp_cell)
lam <- FO$kde$main

BANDS <- list(c(375, 750), c(750, 1500), c(1500, 3000), c(3000, 6000))
NSECT <- 18L                                  # 10-degree sectors over 0-180
NSIM  <- as.integer(Sys.getenv("NSIM_ANISO", 199L))
RMAX  <- max(vapply(BANDS, `[`, numeric(1), 2))

# Axial summaries of a set of orientations, in radians on [0, pi).
#
# TWO harmonics, and the second one matters more than I expected. R2 is the usual
# axial concentration and detects a single preferred axis. It is structurally
# BLIND to four-fold structure: a rose with equal spikes at 0, 45, 90 and 135
# degrees has R2 near zero, because the doubled angles cancel. R4 detects exactly
# that case. The rose diagram showed four-fold spikes at short range that R2
# reported as isotropic, which is why both are computed here.
# THREE harmonics. Each is blind to the next one up, which I learned the hard way
# twice: R2 scored the 375-750 m band as isotropic while the rose showed obvious
# spikes, and R4 then scored the 750-1500 m band as isotropic while the rose still
# showed eight lobes. Spikes at 0, 45, 90 and 135 degrees have a 45-degree period,
# so the 45/135 pair sits in antiphase to the 0/90 pair under exp(4i*theta) and the
# two cancel. exp(8i*theta) is what captures that.
axial <- function(theta) {
  z2 <- mean(exp(2i * theta)); z4 <- mean(exp(4i * theta)); z8 <- mean(exp(8i * theta))
  list(R2 = Mod(z2), axis_deg  = (Arg(z2) / 2) %% pi * 180 / pi,
       R4 = Mod(z4), axis4_deg = (Arg(z4) / 4) %% (pi/2) * 180 / pi,
       R8 = Mod(z8), axis8_deg = (Arg(z8) / 8) %% (pi/4) * 180 / pi)
}

# Orientation of every close pair, folded onto [0, pi).
pair_orient <- function(Y, rmax) {
  cp <- closepairs(Y, rmax = rmax, what = "all")
  keep <- cp$i < cp$j
  list(theta = atan2(cp$dy[keep], cp$dx[keep]) %% pi, d = cp$d[keep])
}

sector_props <- function(theta) {
  b <- pmin(floor(theta / (pi / NSECT)) + 1L, NSECT)
  tabulate(b, nbins = NSECT) / max(length(theta), 1L)
}

message(glue("computing observed orientations, n = {npoints(X)}, rmax = {RMAX} m"))
obs <- pair_orient(X, RMAX)
message(glue("  close pairs: {format(length(obs$theta), big.mark=',')}"))

message(glue("simulating {NSIM} inhomogeneous Poisson patterns"))
sims <- vector("list", NSIM)
for (b in seq_len(NSIM)) {
  Y <- rpoispp(lam)
  sims[[b]] <- pair_orient(Y, RMAX)
  if (b %% 50 == 0) message(glue("  {b}/{NSIM}"))
}

# --- the test, band by band ------------------------------------------------
results <- map_dfr(BANDS, function(bn) {
  r1 <- bn[1]; r2 <- bn[2]
  th_obs <- obs$theta[obs$d >= r1 & obs$d < r2]
  if (length(th_obs) < 50) return(NULL)
  p_obs <- sector_props(th_obs)
  a_obs <- axial(th_obs)

  sim_props <- matrix(NA_real_, nrow = NSIM, ncol = NSECT)
  sim_R <- numeric(NSIM); sim_R4 <- numeric(NSIM); sim_R8 <- numeric(NSIM)
  for (b in seq_len(NSIM)) {
    th <- sims[[b]]$theta[sims[[b]]$d >= r1 & sims[[b]]$d < r2]
    if (!length(th)) { sim_props[b, ] <- rep(1 / NSECT, NSECT); next }
    sim_props[b, ] <- sector_props(th)
    a <- axial(th); sim_R[b] <- a$R2; sim_R4[b] <- a$R4; sim_R8[b] <- a$R8
  }
  null_mean <- colMeans(sim_props)

  # Two statistics. The sector deviation asks whether ANY direction is over- or
  # under-represented; the resultant length asks whether there is a single
  # preferred axis. They can disagree, and if they do that is informative.
  dev_obs <- max(abs(p_obs - null_mean))
  dev_sim <- apply(sim_props, 1, function(p) max(abs(p - null_mean)))
  p_dev <- (1 + sum(dev_sim >= dev_obs)) / (NSIM + 1)
  p_R   <- (1 + sum(sim_R  >= a_obs$R2)) / (NSIM + 1)
  p_R4  <- (1 + sum(sim_R4 >= a_obs$R4)) / (NSIM + 1)
  p_R8  <- (1 + sum(sim_R8 >= a_obs$R8)) / (NSIM + 1)

  tibble(r1 = r1, r2 = r2, n_pairs = length(th_obs),
         R_obs = a_obs$R2, R_null_mean = mean(sim_R),
         R4_obs = a_obs$R4, R4_null_mean = mean(sim_R4),
         R8_obs = a_obs$R8, R8_null_mean = mean(sim_R8),
         axis_deg = a_obs$axis_deg, axis4_deg = a_obs$axis4_deg,
         max_dev = dev_obs, p_sector = p_dev, p_resultant = p_R,
         p_fourfold = p_R4, p_eightfold = p_R8,
         props = list(tibble(sector = seq_len(NSECT),
                             lo_deg = (seq_len(NSECT) - 1) * 180 / NSECT,
                             obs = p_obs, null = null_mean,
                             null_lo = apply(sim_props, 2, quantile, 0.025),
                             null_hi = apply(sim_props, 2, quantile, 0.975))))
})
print(as.data.frame(results |> select(-props) |>
        mutate(across(where(is.numeric), ~ signif(.x, 4)))))

# --- control: is the short-range structure an artefact of my own thinning? --
# ppp_cell keeps at most one point per 375 m cell, so pair displacements at short
# range are dominated by lattice vectors: (1,0) at 375 m, (1,1) at 530 m, (2,0)
# at 750 m. That would put spikes at 0, 45, 90 and 135 degrees whether or not
# fire has any directional structure. ppp_all keeps every detection at its own
# sensor coordinate and has no such lattice, so running the identical test on it
# separates my construction from the data.
message("control: same test on the un-thinned pattern")
XA <- unmark(P$ppp_all)
lam_all <- density(XA, sigma = FO$sigma_main, edge = TRUE, diggle = TRUE, positive = TRUE)
RMAX_C <- 1500; NSIM_C <- 99L
obs_c <- pair_orient(XA, RMAX_C)
th_c  <- obs_c$theta[obs_c$d >= 375 & obs_c$d < 750]
a_c   <- axial(th_c)
simR4_c <- numeric(NSIM_C); simR2_c <- numeric(NSIM_C); simR8_c <- numeric(NSIM_C)
for (b in seq_len(NSIM_C)) {
  Y <- rpoispp(lam_all)
  oc <- pair_orient(Y, RMAX_C)
  tt <- oc$theta[oc$d >= 375 & oc$d < 750]
  if (!length(tt)) next
  a <- axial(tt); simR2_c[b] <- a$R2; simR4_c[b] <- a$R4; simR8_c[b] <- a$R8
}
control <- tibble(
  pattern = c("ppp_cell (one point per 375 m cell)", "ppp_all (every detection)"),
  n = c(npoints(X), npoints(XA)),
  n_pairs = c(sum(obs$d >= 375 & obs$d < 750), length(th_c)),
  R4_obs = c(results$R4_obs[results$r1 == 375], a_c$R4),
  R4_null = c(results$R4_null_mean[results$r1 == 375], mean(simR4_c)),
  p_fourfold = c(results$p_fourfold[results$r1 == 375],
                 (1 + sum(simR4_c >= a_c$R4)) / (NSIM_C + 1)),
  R8_obs = c(results$R8_obs[results$r1 == 375], a_c$R8),
  R8_null = c(results$R8_null_mean[results$r1 == 375], mean(simR8_c)),
  p_eightfold = c(results$p_eightfold[results$r1 == 375],
                  (1 + sum(simR8_c >= a_c$R8)) / (NSIM_C + 1)))
print(as.data.frame(control |> mutate(across(where(is.numeric), ~ signif(.x, 4)))))

# --- does any preferred axis match the canal and road network? -------------
# A7 speculated that fires align with canals. The covariate extract is already
# cached, so the speculation is checkable rather than decorative.
osm <- readRDS(file.path(DIR_DERIVED, "osm_covariates.rds"))
line_orient <- function(lines) {
  if (is.null(lines) || !NROW(lines)) return(NULL)
  g <- suppressWarnings(st_intersection(
    st_geometry(st_transform(lines, STUDY_CRS)), st_geometry(P$win_p)))
  g <- g[!st_is_empty(g)]
  g <- suppressWarnings(st_cast(g, "LINESTRING", warn = FALSE))
  if (!length(g)) return(NULL)
  segs <- do.call(rbind, lapply(seq_along(g), function(k) {
    m <- st_coordinates(g[[k]])[, 1:2, drop = FALSE]
    if (nrow(m) < 2) return(NULL)
    cbind(m[-nrow(m), , drop = FALSE], m[-1, , drop = FALSE])
  }))
  if (is.null(segs)) return(NULL)
  dx <- segs[, 3] - segs[, 1]; dy <- segs[, 4] - segs[, 2]
  len <- sqrt(dx^2 + dy^2)
  th <- atan2(dy, dx) %% pi
  ok <- len > 0
  list(theta = th[ok], len = len[ok])
}
net_summary <- function(o, label) {
  if (is.null(o)) return(NULL)
  # length-weighted, so a long canal counts more than a short kerb
  w <- o$len / sum(o$len)
  z <- sum(w * exp(2i * o$theta))
  b <- pmin(floor(o$theta / (pi / NSECT)) + 1L, NSECT)
  props <- tapply(w, factor(b, levels = seq_len(NSECT)), sum)
  props[is.na(props)] <- 0
  tibble(network = label, R = Mod(z), axis_deg = (Arg(z) / 2) %% pi * 180 / pi,
         total_km = sum(o$len) / 1000,
         props = list(tibble(sector = seq_len(NSECT),
                             lo_deg = (seq_len(NSECT) - 1) * 180 / NSECT,
                             w = as.numeric(props))))
}
networks <- bind_rows(net_summary(line_orient(osm$water), "waterways"),
                      net_summary(line_orient(osm$roads), "roads"))
print(as.data.frame(networks |> select(-props) |>
        mutate(across(where(is.numeric), ~ signif(.x, 4)))))

saveRDS(list(results = results, networks = networks, control = control,
             nsim = NSIM, nsim_control = NSIM_C, nsect = NSECT, bands = BANDS),
        file.path(DIR_DERIVED, "anisotropy.rds"))
message("anisotropy analysis cached")
