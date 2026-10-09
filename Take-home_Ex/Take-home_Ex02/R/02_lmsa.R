# ---------------------------------------------------------------------------
# 02_lmsa.R - spatial weights, global Moran's I, Local Moran's I and Gi*
#
# Output: data/derived/lmsa.rds
#   nb, wt, nb_self, wt_self   the weights used everywhere downstream
#   w_info                     what the weights look like, and the island fix
#   global                     global Moran's I per variable (permutation)
#   lisa                       township x variable Local Moran results
#   gi                         township x variable Gi* results
#   sens                       agreement under raw counts and under k-NN weights
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(sfdep); library(spdep); library(purrr)
})
source("R/00_config.R")
P <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
tsp <- P$tsp
VARS <- c("pv", "battles", "erv", "vac", "fatalities")

## 1. Weights -----------------------------------------------------------------
geom <- st_geometry(tsp)
nb_q <- st_contiguity(geom, queen = TRUE)
islands <- which(card(nb_q) == 0)
# Townships with no shared boundary are given their nearest neighbour by
# centroid distance. Dropping them would remove real places from the map, and
# leaving them with an empty neighbour set gives a zero-weight row that the
# permutation tests cannot use.
cent <- st_centroid(geom)
nb_knn1 <- st_knn(cent, k = K_ISLAND)
nb <- nb_q
for (i in islands) nb[[i]] <- nb_knn1[[i]]
# st_knn links are one-way; make them symmetric so i neighbours j iff j neighbours i.
nb <- spdep::make.sym.nb(nb)
class(nb) <- c("nb", "list")
wt <- st_weights(nb, style = "W")
nb_self <- include_self(nb)
wt_self <- st_weights(nb_self, style = "W")

w_info <- list(
  islands = tibble(TS = tsp$TS[islands], ST = tsp$ST[islands],
                   linked_to = map_chr(islands, \(i) paste(tsp$TS[nb[[i]]], collapse = ", "))),
  card = card(nb), n_links = sum(card(nb)),
  subgraphs = n.comp.nb(nb)$nc
)

## 2. Transformation ---------------------------------------------------------
# Township counts are heavily right-skewed (median 79 political-violence events,
# maximum 1,039; fatalities are dominated by a handful of mass-casualty events).
# Both statistics are cross-products of deviations from the mean, so on raw
# counts a few extreme townships decide the result for their whole
# neighbourhood. log1p() keeps the ordering, keeps zeros as zeros, and lets the
# statistic describe the bulk of townships. Raw counts are kept as a check.
tx <- function(v) log1p(v)

## 3. Global Moran's I ---------------------------------------------------------
set.seed(SEED)
global <- map_dfr(VARS, function(v) {
  g <- global_moran_perm(tx(tsp[[v]]), nb, wt, nsim = NSIM_LMSA)
  tibble(var = v, I = unname(g$statistic), p_sim = g$p.value)
})

## 4. Local Moran's I ---------------------------------------------------------
quad_lab <- c("Low-Low", "High-Low", "Low-High", "High-High")
run_lisa <- function(x, nb, wt) {
  set.seed(SEED)
  lm <- local_moran(x, nb, wt, nsim = NSIM_LMSA)
  tibble(ii = lm$ii, z_ii = lm$z_ii, p_sim = lm$p_ii_sim, quad = as.character(lm$mean))
}
lisa <- map_dfr(VARS, function(v) {
  run_lisa(tx(tsp[[v]]), nb, wt) |>
    mutate(var = v, TS_PCODE = tsp$TS_PCODE, value = tsp[[v]],
           sig = p_sim < ALPHA,
           cls = if_else(sig, quad, "Not significant"),
           p_fdr = p.adjust(p_sim, "fdr"))
})

## 5. Getis-Ord Gi* ------------------------------------------------------------
run_gi <- function(x, nb_self, wt_self) {
  set.seed(SEED)
  g <- local_gstar_perm(x, nb_self, wt_self, nsim = NSIM_LMSA)
  tibble(gi = g$gi_star, p_sim = g$p_sim)
}
gi <- map_dfr(VARS, function(v) {
  run_gi(tx(tsp[[v]]), nb_self, wt_self) |>
    mutate(var = v, TS_PCODE = tsp$TS_PCODE, value = tsp[[v]],
           sig = p_sim < ALPHA,
           cls = case_when(!sig ~ "Not significant",
                           gi > 0 ~ "Hot spot", TRUE ~ "Cold spot"),
           p_fdr = p.adjust(p_sim, "fdr"))
})

## 6. Sensitivity ----------------------------------------------------------------
# (a) raw counts instead of log1p, (b) k-nearest-neighbour weights instead of
# queen contiguity. Agreement is the share of townships given the same class.
nb_k <- st_knn(cent, k = K_SENS); nb_k <- spdep::make.sym.nb(nb_k); class(nb_k) <- c("nb", "list")
wt_k <- st_weights(nb_k)
nb_k_self <- include_self(nb_k); wt_k_self <- st_weights(nb_k_self)

sens <- map_dfr(VARS, function(v) {
  base_l <- filter(lisa, var == v)$cls
  base_g <- filter(gi, var == v)$cls
  raw_l <- run_lisa(tsp[[v]], nb, wt) |> mutate(cls = if_else(p_sim < ALPHA, quad, "Not significant"))
  knn_l <- run_lisa(tx(tsp[[v]]), nb_k, wt_k) |> mutate(cls = if_else(p_sim < ALPHA, quad, "Not significant"))
  raw_g <- run_gi(tsp[[v]], nb_self, wt_self) |>
    mutate(cls = case_when(p_sim >= ALPHA ~ "Not significant", gi > 0 ~ "Hot spot", TRUE ~ "Cold spot"))
  knn_g <- run_gi(tx(tsp[[v]]), nb_k_self, wt_k_self) |>
    mutate(cls = case_when(p_sim >= ALPHA ~ "Not significant", gi > 0 ~ "Hot spot", TRUE ~ "Cold spot"))
  tibble(var = v,
         lisa_raw = mean(base_l == raw_l$cls), lisa_knn = mean(base_l == knn_l$cls),
         gi_raw = mean(base_g == raw_g$cls), gi_knn = mean(base_g == knn_g$cls),
         hh_base = sum(base_l == "High-High"), hh_raw = sum(raw_l$cls == "High-High"),
         hh_knn = sum(knn_l$cls == "High-High"))
})

saveRDS(list(nb = nb, wt = wt, nb_self = nb_self, wt_self = wt_self, w_info = w_info,
             global = global, lisa = lisa, gi = gi, sens = sens),
        file.path(DIR_DERIVED, "lmsa.rds"))

print(w_info$islands); cat("subgraphs:", w_info$subgraphs, " links:", w_info$n_links,
                            " median card:", median(w_info$card), "\n")
print(global)
print(lisa |> count(var, cls) |> pivot_wider(names_from = cls, values_from = n))
print(gi |> count(var, cls) |> pivot_wider(names_from = cls, values_from = n))
print(lisa |> group_by(var) |> summarise(sig = sum(sig), fdr = sum(p_fdr < ALPHA)))
print(sens)
