# ---------------------------------------------------------------------------
# 14 - Residual K for the model that now fits best
#
# R/06 computed a residual K function for its best model, which was distance to
# road plus distance to water. Land cover changed which model that is, so the
# residual has to be recomputed against the new one. Reporting the old residual
# beside the new model would be showing the diagnostic for a model I no longer
# claim.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(purrr)
  library(spatstat.geom); library(spatstat.explore); library(spatstat.model)
})
source("R/00_config.R")
set.seed(SEED)
P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
CM <- readRDS(file.path(DIR_DERIVED, "covariates_models.rds"))
BC <- readRDS(file.path(DIR_DERIVED, "cover_boot.rds"))
X  <- unmark(P$ppp_cell); win <- P$win

lu <- readRDS(file.path(DIR_DERIVED, "landcover.rds")) |>
  st_transform(STUDY_CRS) |> st_make_valid()
lu <- lu[lengths(st_intersects(lu, P$win_p)) > 0, ]
GROUPS <- list(forest = "forest",
               plantation = c("orchard", "farmland", "farmyard", "vineyard"))
dist_to_group <- function(classes) {
  g <- lu |> filter(fclass %in% classes)
  if (!nrow(g)) return(NULL)
  as.im(distfun(as.owin(st_union(st_geometry(g)))), W = win, dimyx = c(600, 300))
}
D_cover <- compact(map(GROUPS, dist_to_group))
dr <- CM$D_road; dw <- CM$D_water
Droad_km <- eval.im(dr / 1000); Dwater_km <- eval.im(dw / 1000)
covs <- list(Droad_km = Droad_km, Dwater_km = Dwater_km)
for (nm in names(D_cover)) { di <- D_cover[[nm]]; covs[[paste0("D", nm, "_km")]] <- eval.im(di / 1000) }

best_name <- BC$aic2$model[1]
form <- switch(best_name,
  all   = "X ~ Droad_km + Dwater_km + Dforest_km + Dplantation_km",
  cover = "X ~ Dforest_km + Dplantation_km",
  access= "X ~ Droad_km + Dwater_km",
  road  = "X ~ Droad_km",
  "X ~ 1")
message("recomputing residual K for the best model: ", best_name, "  (", form, ")")
fit <- ppm(as.formula(form), data = covs)
res_K <- Kres(fit, correction = "border")

saveRDS(list(best_name = best_name, formula = form,
             res_K = as.data.frame(res_K)),
        file.path(DIR_DERIVED, "residual_k.rds"))
message("residual K cached")
