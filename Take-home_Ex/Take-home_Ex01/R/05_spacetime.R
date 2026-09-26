# ---------------------------------------------------------------------------
# 05 - Spatio-temporal structure
#
# A purely spatial analysis of an 87-day fire episode throws away the variable
# that management actually responds to. Three questions are asked here.
#
#  Q1  Is there space-time interaction, i.e. are events that are close in space
#      also close in time by more than the two marginal patterns imply? The null
#      is random labelling: spatial locations are held exactly as observed and
#      only the time labels are permuted. This conditions on both marginals, so
#      a rejection is evidence of genuine interaction rather than of a seasonal
#      trend or a spatial gradient that is already known to exist.
#
#  Q2  Does the spatial signature of the episode change between onset and peak?
#      Answered with a relative-risk surface, which is a ratio of two kernel
#      estimates and so is interpretable as "where is the peak concentrated
#      relative to the onset", not merely "where are there more fires".
#
#  Q3  Does the episode move? Answered with the weekly mean centre and standard
#      distance, which is a crude but honest summary that a duty officer could
#      read off a map.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(spatstat.geom); library(spatstat.explore)
  library(spatstat.random)
})
source("R/00_config.R")
set.seed(SEED)
P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
FO <- readRDS(file.path(DIR_DERIVED, "firstorder.rds"))

XD <- P$ppp_daycell     # the spatio-temporal unit: one event per cell per day
tt <- as.integer(marks(XD)$doy)
message(sprintf("ppp_daycell n = %d over %d days", npoints(XD), P$n_days))

# --- Q1. Knox test of space-time interaction -------------------------------
# closepairs is computed once at the largest spatial threshold; the permutation
# then only reshuffles time labels, so 999 permutations are affordable even at
# this sample size. Reported across a grid of (s, t) because a single pair of
# thresholds chosen after seeing the answer would not be credible.
S_THRESH <- c(500, 1000, 2000, 5000)      # metres
T_THRESH <- c(1, 3, 7, 14)                # days

# Pairs are enumerated once. Only the time labels are permuted, so the spatial
# part of the computation is never repeated.
cp <- closepairs(XD, rmax = max(S_THRESH), what = "all")
keep <- cp$i < cp$j                        # each unordered pair exactly once
i <- cp$i[keep]; j <- cp$j[keep]; d_sp <- cp$d[keep]
rm(cp); gc(verbose = FALSE)
message(sprintf("close pairs within %d m: %s", max(S_THRESH),
                format(length(i), big.mark = ",")))

# The naive loop would evaluate every (s, t) cell for every permutation, which at
# this pair count is tens of billions of comparisons. Instead each pair is binned
# once in space and once per permutation in time, the joint 2-D histogram is
# tabulated, and a double cumulative sum recovers all (s, t) cells at once. The
# result is identical; the cost falls by two orders of magnitude.
nS <- length(S_THRESH); nT <- length(T_THRESH)
# s_bin = index of the smallest threshold the pair satisfies (nS+1 = satisfies none)
s_bin <- findInterval(d_sp, S_THRESH, left.open = TRUE) + 1L
rm(d_sp); gc(verbose = FALSE)

knox_counts <- function(tvec) {
  dt <- abs(tvec[i] - tvec[j])
  t_bin <- findInterval(dt, T_THRESH, left.open = TRUE) + 1L
  tab <- tabulate(s_bin + (nS + 1L) * (t_bin - 1L), nbins = (nS + 1L) * (nT + 1L))
  m <- matrix(tab, nrow = nS + 1L, ncol = nT + 1L)
  # cumulative in both directions: cell [a, b] = pairs with s_bin<=a and t_bin<=b
  cm <- apply(apply(m, 2, cumsum), 1, cumsum)   # returns t x s
  t(cm)[seq_len(nS), seq_len(nT), drop = FALSE]
}

NPERM    <- as.integer(Sys.getenv("NPERM_KNOX", 999L))
obs_mat  <- knox_counts(tt)
sim_arr  <- array(0L, dim = c(nS, nT, NPERM))
for (b in seq_len(NPERM)) sim_arr[, , b] <- knox_counts(sample(tt))

knox_tab <- expand.grid(si = seq_len(nS), ti = seq_len(nT)) |>
  as_tibble() |>
  mutate(s_m = S_THRESH[si], t_days = T_THRESH[ti],
         observed = mapply(function(a, b) obs_mat[a, b], si, ti),
         exp_perm = mapply(function(a, b) mean(sim_arr[a, b, ]), si, ti),
         ratio = observed / exp_perm,
         # Monte Carlo p-value including the observed value, so p is never zero.
         p_value = mapply(function(a, b)
           (1 + sum(sim_arr[a, b, ] >= obs_mat[a, b])) / (NPERM + 1), si, ti)) |>
  select(s_m, t_days, observed, exp_perm, ratio, p_value) |>
  arrange(s_m, t_days)
print(as.data.frame(knox_tab))
n_close_pairs <- length(i)

# --- Q2. relative risk, onset versus peak ----------------------------------
# Split at the change point in the daily series rather than at a round date.
# The change point is taken as the first day on which the 7-day mean exceeds
# half its eventual maximum, which is a rule stated before the surfaces are
# inspected.
roll <- FO$daily$roll7
cp_day <- which(roll >= 0.5 * max(roll, na.rm = TRUE))[1]
split_date <- WINDOW_START + cp_day - 1L
message(sprintf("onset/peak split at day %d (%s)", cp_day, split_date))

season <- factor(ifelse(marks(XD)$date_local < split_date, "onset", "peak"),
                 levels = c("onset", "peak"))
XS <- XD; marks(XS) <- season
XS <- XS[!is.na(marks(XS))]

# bw.relrisk chooses a common bandwidth for the two components by cross-
# validation of the risk surface itself, which is the right target here; using
# each type's own optimal bandwidth would make the ratio uninterpretable.
sig_rr <- tryCatch(as.numeric(bw.relrisk(XS, method = "likelihood")),
                   error = function(e) FO$sigma_main)
rr <- relrisk(XS, sigma = sig_rr, casecontrol = TRUE, case = "peak",
              relative = FALSE, edge = TRUE)
# Pointwise tolerance contours: where does the probability of belonging to the
# peak period depart from the overall share by more than sampling error?
rr_se <- relrisk(XS, sigma = sig_rr, casecontrol = TRUE, case = "peak",
                 relative = FALSE, edge = TRUE, se = TRUE)

# --- Q3. weekly mean centre and standard distance -------------------------
co <- coords(XD)
traj <- tibble(x = co$x, y = co$y, doy = tt,
               week = ((tt - 1) %/% 7) + 1L) |>
  group_by(week) |>
  summarise(n = n(), mx = mean(x), my = mean(y),
            sd_dist = sqrt(mean((x - mean(x))^2 + (y - mean(y))^2)),
            .groups = "drop") |>
  mutate(week_start = WINDOW_START + (week - 1) * 7,
         step_km = c(NA, sqrt(diff(mx)^2 + diff(my)^2) / 1000),
         sd_km = sd_dist / 1000)
print(as.data.frame(traj |> select(week, week_start, n, sd_km, step_km)))

# --- Q4. separable space-time intensity, by weighted KDE -------------------
# A product-kernel estimate: a fixed spatial bandwidth with a Gaussian weight in
# time. Implemented directly rather than via a dedicated package, so the
# separability assumption is explicit: the spatial kernel does not change with
# time, only the weight each event receives.
st_slices <- seq(7, P$n_days - 7, by = 14)
sigma_t <- 5   # days
st_kde <- lapply(st_slices, function(tc) {
  w <- dnorm(tt - tc, sd = sigma_t)
  im <- density(XD, sigma = FO$sigma_main, weights = w / sum(w),
                edge = TRUE, diggle = TRUE, positive = TRUE)
  list(day = tc, date = WINDOW_START + tc - 1L, im = im,
       n_eff = sum(w)^2 / sum(w^2))
})

saveRDS(list(
  knox_tab = knox_tab, n_close_pairs = n_close_pairs, nperm = NPERM,
  S_THRESH = S_THRESH, T_THRESH = T_THRESH,
  split_date = split_date, cp_day = cp_day, sig_rr = sig_rr,
  rr = rr, rr_se = rr_se, season_counts = table(season),
  traj = traj, st_kde = st_kde, st_slices = st_slices, sigma_t = sigma_t
), file.path(DIR_DERIVED, "spacetime.rds"))
message("spatio-temporal analysis cached")
