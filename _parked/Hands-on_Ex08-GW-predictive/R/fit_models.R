# ---------------------------------------------------------------------------
# fit_models.R - every slow step of Hands-on Exercise 8, cached to data/model/
#
# The page shows these chunks with eval: false and reads the cached objects,
# the same pattern the handout uses. Each step is skipped if its file exists.
# Run from Hands-on_Ex/Hands-on_Ex08/:  Rscript R/fit_models.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(GWmodel); library(SpatialML); library(ranger)
  library(rsample); library(tidyverse)
})
dir.create("data/model", showWarnings = FALSE)
p <- function(f) file.path("data/model", f)
tick <- function(label, expr) {
  t0 <- Sys.time(); force(expr)
  cat(sprintf("%-28s %6.1f min\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

mdata <- read_rds("data/rawdata/mdata.rds")
if (!file.exists(p("train_data.rds"))) {
  set.seed(1234)
  resale_split <- initial_split(mdata, prop = 6.5/10,)
  write_rds(training(resale_split), p("train_data.rds"))
  write_rds(testing(resale_split), p("test_data.rds"))
}
train_data <- read_rds(p("train_data.rds")); test_data <- read_rds(p("test_data.rds"))

f <- resale_price ~ floor_area_sqm + storey_order + remaining_lease_mths +
  PROX_CBD + PROX_ELDERLYCARE + PROX_HAWKER + PROX_MRT + PROX_PARK + PROX_MALL +
  PROX_SUPERMARKET + WITHIN_350M_KINDERGARTEN + WITHIN_350M_CHILDCARE +
  WITHIN_350M_BUS + WITHIN_1KM_PRISCH

if (!file.exists(p("price_mlr.rds"))) write_rds(lm(f, data = train_data), p("price_mlr.rds"))

coords_train <- st_coordinates(train_data); coords_test <- st_coordinates(test_data)
write_rds(coords_train, p("coords_train.rds")); write_rds(coords_test, p("coords_test.rds"))
train_nogeom <- st_drop_geometry(train_data)

if (!file.exists(p("rf.rds"))) tick("ranger", {
  set.seed(1234); write_rds(ranger(f, data = train_nogeom), p("rf.rds"))
})

if (!file.exists(p("bw_adaptive.rds"))) tick("bw.gwr (train)", {
  bw <- bw.gwr(f, data = train_data, approach = "CV", kernel = "gaussian",
               adaptive = TRUE, longlat = FALSE)
  write_rds(bw, p("bw_adaptive.rds"))
})
bw_adaptive <- read_rds(p("bw_adaptive.rds")); cat("bw_adaptive =", bw_adaptive, "\n")

if (!file.exists(p("gwr_adaptive.rds"))) tick("gwr.basic", {
  write_rds(gwr.basic(f, data = train_data, bw = bw_adaptive, kernel = "gaussian",
                      adaptive = TRUE, longlat = FALSE), p("gwr_adaptive.rds"))
})

if (!file.exists(p("gwr_pred.rds"))) tick("GWR prediction", {
  # The handout uses gwr.predict(). In GWmodel 2.4.1 that call fails for more
  # than 5,000 training points (a positional-argument bug in its distance
  # loop), and supplying the distance matrices to get round it leaves it
  # computing a prediction variance through t(I - S) %*% (I - S) on
  # 10,335 x 10,335 dense matrices, which did not finish in over an hour.
  # The predictions themselves are only local coefficients at each test
  # location times that flat's attributes, so they are computed directly:
  # gwr.basic() calibrates on the training flats at the test locations.
  # R/check_gwr_predict.R shows this matches gwr.predict() on a subset small
  # enough for gwr.predict() to finish.
  g <- gwr.basic(f, data = train_data, regression.points = test_data, bw = 40,
                 kernel = "gaussian", adaptive = TRUE, longlat = FALSE)
  B <- as.matrix(st_drop_geometry(st_as_sf(g$SDF))[, 1:15])
  X <- cbind(1, as.matrix(st_drop_geometry(test_data)[, all.vars(f)[-1]]))
  write_rds(tibble(prediction = rowSums(B * X)), p("gwr_pred.rds"))
})

if (!file.exists(p("gwRF_adaptive.rds"))) tick("grf", {
  set.seed(1234)
  g <- grf(f, dframe = train_nogeom, bw = 40, kernel = "adaptive", coords = coords_train)
  write_rds(g, p("gwRF_adaptive.rds"))
})

if (!file.exists(p("GRF_pred.rds"))) tick("predict.grf", {
  g <- read_rds(p("gwRF_adaptive.rds"))
  test_nogeom <- cbind(test_data, coords_test) |> st_drop_geometry()
  # SpatialML 1.8.2 no longer exports predict.grf(); it is an S3 method,
  # so the generic predict() dispatches to it.
  write_rds(predict(g, test_nogeom, x.var.name = "X", y.var.name = "Y",
                        local.w = 1, global.w = 0), p("GRF_pred.rds"))
})

# The handout also computes a bandwidth on the test data. It is never used
# (prediction is fixed at bw = 40), and choosing a bandwidth on the test set
# would let the test data shape the model it judges, so it is not computed.
cat("done\n")
