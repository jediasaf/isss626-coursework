# ---------------------------------------------------------------------------
# spatial_holdout.R - how much of the test accuracy is location leakage?
#
# The handout's split is random over flats. Many flats share a block, so most
# test flats have a training flat at the same coordinates. This script repeats
# the comparison with a spatially blocked holdout: whole 2 km grid cells go to
# the test set, so no test flat has a training neighbour in the same block.
# Five different random draws of held-out cells, for MLR, random forest and
# GWR (bandwidth kept at the 40 neighbours chosen on the random split, and
# predicted with gwr.basic(regression.points = ...) as in fit_models.R). The
# geographically weighted random forest is left out: at about nine minutes per
# fit and prediction it would add three quarters of an hour.
# Output: data/model/holdout.rds
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(sf); library(ranger); library(GWmodel); library(tidyverse); library(Metrics)})
tr <- read_rds("data/model/train_data.rds"); te <- read_rds("data/model/test_data.rds")
f <- formula(read_rds("data/model/price_mlr.rds"))
all <- bind_rows(tr, te)
xy <- st_coordinates(all)
cell <- paste(floor(xy[, 1] / 2000), floor(xy[, 2] / 2000))
cells <- unique(cell)

nn_dist <- function(a, b) as.numeric(st_distance(a, b[st_nearest_feature(a, b), ], by_element = TRUE))

blocked <- map_dfr(1:5, function(s) {
  set.seed(s)
  hold <- sample(cells, round(0.35 * length(cells)))
  btr <- all[!cell %in% hold, ]; bte <- all[cell %in% hold, ]
  m <- lm(f, btr)
  set.seed(1234); r <- ranger(f, st_drop_geometry(btr))
  g <- gwr.basic(f, data = btr, regression.points = bte, bw = 40, kernel = "gaussian",
                 adaptive = TRUE, longlat = FALSE)
  B <- as.matrix(st_drop_geometry(st_as_sf(g$SDF))[, 1:15])
  X <- cbind(1, as.matrix(st_drop_geometry(bte)[, all.vars(f)[-1]]))
  tibble(draw = s, n_train = nrow(btr), n_test = nrow(bte),
         median_nn_m = median(nn_dist(bte, btr)),
         MLR = rmse(bte$resale_price, predict(m, bte)),
         RF = rmse(bte$resale_price, predict(r, st_drop_geometry(bte))$predictions),
         GWR = rmse(bte$resale_price, rowSums(B * X)))
})

leak <- tibble(
  same_location = mean(nn_dist(te, tr) < 1),
  within_100m = mean(nn_dist(te, tr) < 100),
  n_cells = length(cells))

write_rds(list(blocked = blocked, leak = leak), "data/model/holdout.rds")
print(leak); print(blocked)
