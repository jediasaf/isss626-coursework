# ---------------------------------------------------------------------------
# 15 - Spatial autocorrelation on an areal aggregation, and why the obvious
#      version of it is uninformative
#
# Everything so far treats the data as a point pattern. Aggregating to areal
# units and running Moran's I is a different lens on the same question, and it is
# the lens the spatial-weights material in this course provides. It is worth
# doing because it is independent machinery: contiguity weights and a
# permutation test on counts share no assumptions with Kinhom and its envelopes.
#
# It has to be set up carefully, because the naive version repeats the mistake
# Section 7.1 is built around. Moran's I on raw burned-cell counts will be large
# and significant here whatever the fire does, because intensity varies by two
# orders of magnitude between the peat belt and the northern interior. Adjacent
# hexagons in the south are both high and adjacent hexagons in the north are both
# low, so counts are autocorrelated by construction. That is the first-order
# trend, measured again.
#
# The informative version asks whether autocorrelation survives the trend. Each
# hexagon gets an expected count by integrating lambda-hat over it, and the
# variable analysed is the Pearson residual
#
#     (observed - expected) / sqrt(expected)
#
# which is the areal analogue of replacing CSR with an inhomogeneous Poisson null.
# Both are reported, because the contrast is the point.
#
# Hexagons, not squares. Section 7.5 found that a square 375 m lattice imprints
# its own geometry on short-range structure, and a square aggregation grid has the
# same defect: four edge neighbours and four corner neighbours at different
# distances. A hexagonal tessellation gives six neighbours all equidistant.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(purrr); library(glue); library(tidyr)
  library(sfdep); library(spdep); library(spatstat.geom)
})
source("R/00_config.R")
set.seed(SEED)

P  <- readRDS(file.path(DIR_DERIVED, "prepared.rds"))
FO <- readRDS(file.path(DIR_DERIVED, "firstorder.rds"))
X  <- P$ppp_cell
win_p <- P$win_p
lam <- FO$kde$main

pts <- st_as_sf(data.frame(x = X$x, y = X$y), coords = c("x", "y"), crs = STUDY_CRS)

# lambda-hat as pixel centres carrying expected counts, so integrating over any
# polygon is a spatial join and a sum instead of 1,200 separate integrals.
lam_df <- as.data.frame(lam)
names(lam_df) <- c("x", "y", "v")
px_area <- lam$xstep * lam$ystep
lam_pts <- st_as_sf(lam_df[is.finite(lam_df$v), ], coords = c("x", "y"), crs = STUDY_CRS) |>
  mutate(expected = v * px_area)

build_hex <- function(cellsize_m) {
  g <- st_make_grid(win_p, cellsize = cellsize_m, square = FALSE) |> st_as_sf()
  names(g)[1] <- "geometry"; st_geometry(g) <- "geometry"
  g <- g[lengths(st_intersects(g, win_p)) > 0, ]
  # keep hexagons that are mostly inside the regency, so edge units are not
  # half-empty and artificially low
  g$frac_in <- as.numeric(st_area(st_intersection(g, win_p))) / as.numeric(st_area(g))
  g <- g[g$frac_in > 0.5, ]
  g$hex_id <- seq_len(nrow(g))
  g$n <- lengths(st_intersects(g, pts))
  e <- st_join(lam_pts, g["hex_id"], join = st_within) |> st_drop_geometry() |>
    filter(!is.na(hex_id)) |> group_by(hex_id) |> summarise(expected = sum(expected), .groups = "drop")
  g <- left_join(g, e, by = "hex_id") |> mutate(expected = replace_na(expected, 0))
  # rescale so total expected equals total observed; the kernel estimate integrates
  # to slightly less than n after edge correction and cropping
  g$expected <- g$expected * sum(g$n) / sum(g$expected)
  g$resid <- (g$n - g$expected) / sqrt(pmax(g$expected, 0.5))
  g
}

run_moran <- function(g, var, label, nsim = 999) {
  nb <- st_contiguity(g, queen = TRUE)
  ok <- lengths(nb) > 0
  if (any(!ok)) { g <- g[ok, ]; nb <- st_contiguity(g, queen = TRUE) }
  wt <- st_weights(nb, style = "W")
  x <- g[[var]]
  gm <- global_moran_test(x, nb, wt)
  list(g = g, nb = nb, wt = wt,
       row = tibble(variable = label, n_units = nrow(g),
                    mean_nb = mean(lengths(nb)),
                    moran_I = unname(gm$estimate[1]),
                    expectation = unname(gm$estimate[2]),
                    z = unname(gm$statistic), p = gm$p.value))
}

# --- global Moran's I across three resolutions, both variables -------------
SIZES <- c(2000, 3000, 5000)
hexes <- set_names(map(SIZES, build_hex), paste0(SIZES / 1000, " km"))
global <- imap_dfr(hexes, function(g, sz) {
  bind_rows(run_moran(g, "n", "raw burned-cell count")$row,
            run_moran(g, "resid", "Pearson residual after lambda-hat")$row) |>
    mutate(hex_size = sz, .before = 1)
})
print(as.data.frame(global |> mutate(across(where(is.numeric), ~ signif(.x, 4)))))

# --- LISA at the middle resolution -----------------------------------------
G <- hexes[["3 km"]]
mr <- run_moran(G, "resid", "resid")
G <- mr$g
lm_raw   <- local_moran(G$n,     mr$nb, mr$wt, nsim = 999)
lm_resid <- local_moran(G$resid, mr$nb, mr$wt, nsim = 999)

lag_of <- function(x) as.numeric(spdep::lag.listw(spdep::nb2listw(mr$nb, style = "W"), x))

classify <- function(x, lm, alpha = 0.05) {
  p <- lm[["p_folded_sim"]]
  # Benjamini-Hochberg. With ~1,200 units an unadjusted 5% test expects about 60
  # false positives, which would paint a convincing but meaningless map.
  padj <- p.adjust(p, method = "BH")
  xs <- as.numeric(scale(x)); ls <- as.numeric(scale(lag_of(x)))
  cl <- case_when(padj >= alpha ~ "not significant",
                  xs > 0 & ls > 0 ~ "High-High",
                  xs < 0 & ls < 0 ~ "Low-Low",
                  xs > 0 & ls < 0 ~ "High-Low",
                  TRUE ~ "Low-High")
  tibble(Ii = lm[["ii"]], p = p, p_adj = padj, cluster = cl)
}
cl_raw   <- classify(G$n, lm_raw)
cl_resid <- classify(G$resid, lm_resid)
G$cl_raw <- cl_raw$cluster; G$cl_resid <- cl_resid$cluster
G$Ii_resid <- cl_resid$Ii; G$padj_resid <- cl_resid$p_adj

counts <- bind_rows(
  count(cl_raw, cluster) |> mutate(variable = "raw count"),
  count(cl_resid, cluster) |> mutate(variable = "residual")) |>
  pivot_wider(names_from = variable, values_from = n, values_fill = 0)
print(as.data.frame(counts))

saveRDS(list(global = global, hex = G, sizes = SIZES,
             counts = counts, alpha = 0.05, nsim = 999,
             n_units = nrow(G), mean_nb = mean(lengths(mr$nb))),
        file.path(DIR_DERIVED, "lisa.rds"))
message("LISA cached")
