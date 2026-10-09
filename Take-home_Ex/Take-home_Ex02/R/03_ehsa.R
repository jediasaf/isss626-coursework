# ---------------------------------------------------------------------------
# 03_ehsa.R - space-time cube, Gi* per quarter, Mann-Kendall, EHSA
#
# Output: data/derived/ehsa.rds
#   gi_q      township x quarter Gi* and permutation p (spatial weights only)
#   mk        township x variable Mann-Kendall tau and p on the Gi* series
#   ehsa      township x variable EHSA classification (sfdep, one-quarter lag)
#   check     the location-order check on sfdep's output
#   sens      classification agreement: raw counts, and a monthly cube
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(sfdep); library(purrr); library(Kendall)
})
source("R/00_config.R")
P <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
L <- readRDS(file.path(DIR_DERIVED, "lmsa.rds"))
tsp <- P$tsp
stopifnot(!is.unsorted(tsp$TS_PCODE))
VARS <- c("pv", "battles", "erv", "vac")

## 1. The cube -------------------------------------------------------------------
# sfdep wants an integer-like time index. Quarters are numbered 1..19 and the
# date of each is kept alongside for labelling.
geo <- tsp |> select(TS_PCODE)
geo$nb <- L$nb_self; geo$wt <- L$wt_self
mk_cube <- function(cube, periods) {
  cube |> mutate(t = match(period, periods)) |> arrange(t, TS_PCODE)
}
cube <- mk_cube(P$cube, P$quarters)
st_q <- spacetime(cube, geo, .loc_col = "TS_PCODE", .time_col = "t")
stopifnot(is_spacetime_cube(st_q))

## 2. Gi* per quarter, then Mann-Kendall --------------------------------------
# This is the brief's "Mann-Kendall test using the spatio-temporal local Gi*":
# compute Gi* for every township in every quarter with the spatial weights,
# then test each township's 19-value Gi* series for a monotonic trend.
tx <- function(v) log1p(v)
set.seed(SEED)
gi_q <- map_dfr(VARS, function(v) {
  cube |> group_by(t) |>
    group_modify(function(d, k) {
      g <- local_gstar_perm(tx(d[[v]]), L$nb_self, L$wt_self, nsim = NSIM_EHSA)
      tibble(TS_PCODE = d$TS_PCODE, value = d[[v]], gi = g$gi_star, p_sim = g$p_sim)
    }) |> ungroup() |> mutate(var = v)
})
gi_q$period <- P$quarters[gi_q$t]

mk <- gi_q |> group_by(var, TS_PCODE) |>
  summarise(mk = list(unclass(MannKendall(gi))), mean_gi = mean(gi),
            n_hot = sum(p_sim < ALPHA & gi > 0), n_cold = sum(p_sim < ALPHA & gi < 0),
            .groups = "drop") |>
  unnest_wider(mk) |>
  mutate(sl = as.numeric(sl), tau = as.numeric(tau),
         trend = case_when(sl >= ALPHA ~ "No significant trend",
                           tau > 0 ~ "Rising", TRUE ~ "Falling"))

## 3. EHSA ---------------------------------------------------------------------
# sfdep's implementation adds a one-period time lag (k = 1): a township's
# neighbourhood in quarter t includes its neighbours in quarter t-1. threshold
# sets the significance used inside the classification rules; it is set to the
# brief's 0.05 rather than sfdep's default of 0.01.
run_ehsa <- function(stc, v, f = tx, nsim = NSIM_EHSA) {
  set.seed(SEED)
  stc[[".x"]] <- f(stc[[v]])
  e <- emerging_hotspot_analysis(stc, ".x", k = 1, nb_col = "nb", wt_col = "wt",
                                 nsim = nsim, threshold = ALPHA)
  # Location-order guard. sfdep 0.2.5 labels its rows with the geometry's
  # location column but computes them in alphabetical order of location. The
  # geometry was sorted by PCode in 01_prepare.R so the two coincide; this
  # assertion makes the script fail rather than mislabel if that ever changes.
  stopifnot(identical(as.character(e$location), sort(unique(as.character(stc$TS_PCODE)))))
  tibble(TS_PCODE = as.character(e$location), tau = e$tau, p_value = e$p_value,
         classification = e$classification, var = v)
}
ehsa <- map_dfr(VARS, \(v) run_ehsa(st_q, v))

## 4. Independent check on the labels ------------------------------------------
# Join EHSA back to the separately computed mean Gi*. If labels are attached to
# the right townships, hot spot classes must have positive mean Gi* and cold
# spot classes negative.
check <- ehsa |> left_join(select(mk, var, TS_PCODE, mean_gi), by = c("var", "TS_PCODE")) |>
  mutate(kind = case_when(grepl("hotspot", classification) ~ "hot",
                          grepl("coldspot", classification) ~ "cold", TRUE ~ "none")) |>
  group_by(var, kind) |> summarise(n = n(), mean_gi = mean(mean_gi), .groups = "drop")
stopifnot(all(check$mean_gi[check$kind == "hot"] > 0),
          all(check$mean_gi[check$kind == "cold"] < 0))

## 5. Sensitivity -----------------------------------------------------------------
# (a) raw counts in place of log1p; (b) a monthly cube (57 periods) in place of
# the quarterly one. Agreement is the share of townships with the same label.
ehsa_raw <- run_ehsa(st_q, "pv", f = identity)
cube_m <- mk_cube(P$cube_m, P$months)
st_m <- spacetime(cube_m, geo, .loc_col = "TS_PCODE", .time_col = "t")
ehsa_m <- run_ehsa(st_m, "pv")
base <- filter(ehsa, var == "pv")
agree <- function(a, b) mean(a$classification == b$classification[match(a$TS_PCODE, b$TS_PCODE)])
sens <- list(
  agree_raw = agree(base, ehsa_raw), agree_month = agree(base, ehsa_m),
  xtab_month = table(quarterly = base$classification,
                     monthly = ehsa_m$classification[match(base$TS_PCODE, ehsa_m$TS_PCODE)]),
  ehsa_m = ehsa_m, ehsa_raw = ehsa_raw)

saveRDS(list(gi_q = gi_q, mk = mk, ehsa = ehsa, check = check, sens = sens),
        file.path(DIR_DERIVED, "ehsa.rds"))

print(ehsa |> count(var, classification) |> pivot_wider(names_from = var, values_from = n), n = 30)
print(mk |> count(var, trend) |> pivot_wider(names_from = trend, values_from = n))
print(check)
cat("agreement raw:", sens$agree_raw, " monthly:", sens$agree_month, "\n")
