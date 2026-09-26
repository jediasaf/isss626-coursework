# ---------------------------------------------------------------------------
# 16 - Which land-cover term actually carries the signal
#
# The covariate section reported forest and plantation together and read the
# coefficient signs as an agricultural frontier. Checking the support behind each
# curve undercuts half of that. The median burned cell sits 19.6 km from mapped
# forest and 4.1 km from mapped plantation, and the forest rhohat curve is not
# monotone. A positive linear coefficient over a range that wide is a crude
# summary of a shape that is not linear.
#
# This fits the two terms separately so the claim can be attached to whichever
# one supports it.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(sf); library(purrr); library(glue)
  library(spatstat.geom); library(spatstat.explore); library(spatstat.model)
})
source("R/00_config.R"); set.seed(SEED)
P <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
CM <- readRDS(file.path(DIR_DERIVED, "covariates_models.rds"))
X <- unmark(P$ppp_cell); win <- P$win

lu <- readRDS(file.path(DIR_DERIVED, "landcover.rds")) |>
  st_transform(STUDY_CRS) |> st_make_valid()
lu <- lu[lengths(st_intersects(lu, P$win_p)) > 0, ]
G <- list(forest = "forest", plantation = c("orchard","farmland","farmyard","vineyard"))
D <- imap(G, function(cl, nm)
  as.im(distfun(as.owin(st_union(st_geometry(lu |> filter(fclass %in% cl))))),
        W = win, dimyx = c(600, 300)))

# how far the data actually sits from each feature
support <- imap_dfr(D, function(d, nm) {
  v <- d[X]
  tibble(feature = nm, median_km = median(v)/1000, p75_km = quantile(v,.75)/1000,
         p95_km = quantile(v,.95)/1000, share_beyond_20km = mean(v > 20000))
})
print(as.data.frame(support |> mutate(across(where(is.numeric), ~ round(.x, 2)))))

dr <- CM$D_road; dw <- CM$D_water
covs <- list(Droad_km = eval.im(dr/1000), Dwater_km = eval.im(dw/1000))
df <- D$forest; dp <- D$plantation
covs$Dforest_km <- eval.im(df/1000); covs$Dplantation_km <- eval.im(dp/1000)

fitm <- function(f) ppm(as.formula(f), data = covs)
mods <- list(
  null        = ppm(X ~ 1),
  road        = fitm("X ~ Droad_km"),
  forest      = fitm("X ~ Dforest_km"),
  plantation  = fitm("X ~ Dplantation_km"),
  cover_both  = fitm("X ~ Dforest_km + Dplantation_km"),
  all         = fitm("X ~ Droad_km + Dwater_km + Dforest_km + Dplantation_km"))
aic <- tibble(model = names(mods), AIC = sapply(mods, AIC),
              npar = sapply(mods, function(m) length(coef(m)))) |>
  mutate(dAIC = AIC - min(AIC), gain_over_null = AIC[model == "null"] - AIC) |>
  arrange(AIC)
print(as.data.frame(aic |> mutate(across(where(is.numeric), ~ round(.x, 1)))))

saveRDS(list(support = support, aic = aic,
             coef_all = coef(summary(mods$all))),
        file.path(DIR_DERIVED, "cover_split.rds"))
message("cover split cached")
