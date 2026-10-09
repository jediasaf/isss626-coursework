# ---------------------------------------------------------------------------
# 03b_ehsa_gi.R - keep the space-time Gi* series EHSA classified on
#
# EHSA classifies each township from its own space-time Gi* (neighbours in the
# current and previous quarter), not from the purely spatial Gi* of 03_ehsa.R.
# To explain a label, that series is needed. This reruns the all-violence EHSA
# with include_gi = TRUE under the same seed and confirms the labels match.
# Output: data/derived/ehsa_gi.rds
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(sf); library(dplyr); library(sfdep)})
source("R/00_config.R")
P <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
L <- readRDS(file.path(DIR_DERIVED, "lmsa.rds"))
E <- readRDS(file.path(DIR_DERIVED, "ehsa.rds"))
geo <- P$tsp |> select(TS_PCODE); geo$nb <- L$nb_self; geo$wt <- L$wt_self
cube <- P$cube |> mutate(t = match(period, P$quarters)) |> arrange(t, TS_PCODE)
stc <- spacetime(cube, geo, .loc_col = "TS_PCODE", .time_col = "t")
set.seed(SEED)
stc[[".x"]] <- log1p(stc[["pv"]])
e <- emerging_hotspot_analysis(stc, ".x", k = 1, nb_col = "nb", wt_col = "wt",
                               nsim = NSIM_EHSA, threshold = ALPHA, include_gi = TRUE)
base <- filter(E$ehsa, var == "pv")
stopifnot(identical(as.character(e$location), base$TS_PCODE),
          identical(e$classification, base$classification))
# gi rows follow the cube order (time, then PCode)
gis <- attr(e, "gi_star")
gi_spt <- tibble(TS_PCODE = cube$TS_PCODE, t = cube$t, period = cube$period,
                 gi = gis$gi_star, p_sim = gis$p_sim)
saveRDS(gi_spt, file.path(DIR_DERIVED, "ehsa_gi.rds"))
osc <- base$TS_PCODE[base$classification == "oscilating hotspot"]
g <- filter(gi_spt, TS_PCODE %in% osc, p_sim < ALPHA, gi < 0)
cat("oscillating:", length(osc), " with sig cold quarter:", n_distinct(g$TS_PCODE), "\n"); print(table(g$t))
