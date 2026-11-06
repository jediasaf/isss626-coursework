# ---------------------------------------------------------------------------
# handout_changes.R - the two steps the chapter 14 handout added or changed on
# 19 October 2026, measured against the cached split in data/model/
#
#   1. overlapping points and st_jitter(amount = 1)
#   2. GWR test predictions at bw = 200 (the handout's gwr.predict() value)
#      as well as bw = 40 (the cross-validated bandwidth)
#
# Run from Hands-on_Ex/Hands-on_Ex10/:  Rscript R/handout_changes.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(GWmodel); library(Metrics); library(tidyverse)
})
p <- function(f) file.path("data/model", f)

mdata <- read_rds("data/rawdata/mdata.rds")
train_data <- read_rds(p("train_data.rds"))
test_data <- read_rds(p("test_data.rds"))

# 1. Overlapping points -------------------------------------------------------
n_overlap <- function(x) sum(lengths(st_equals(x, x)) > 1)

set.seed(1234)
mdata_jit <- st_jitter(mdata, amount = 1)
shift <- as.numeric(st_distance(mdata, mdata_jit, by_element = TRUE))

overlap <- tibble(
  flats = nrow(mdata),
  locations = nrow(distinct(as_tibble(st_coordinates(mdata)))),
  overlapping_before = n_overlap(mdata),
  overlapping_after = n_overlap(mdata_jit),
  max_shift_m = max(shift)
)

# Does jitter remove the twin problem of section 11? Jitter both sets the same
# way and measure each test flat's distance to its nearest training flat.
set.seed(1234)
train_jit <- st_jitter(train_data, amount = 1)
test_jit <- st_jitter(test_data, amount = 1)
nn <- st_nearest_feature(test_jit, train_jit)
d <- as.numeric(st_distance(test_jit, train_jit[nn, ], by_element = TRUE))
twins <- tibble(
  exact_same_location = mean(d == 0),
  within_3m = mean(d <= 3),
  within_100m = mean(d <= 100)
)

# 2. GWR test predictions at the handout's bandwidth ---------------------------
f <- resale_price ~ floor_area_sqm + storey_order + remaining_lease_mths +
  PROX_CBD + PROX_ELDERLYCARE + PROX_HAWKER + PROX_MRT + PROX_PARK + PROX_MALL +
  PROX_SUPERMARKET + WITHIN_350M_KINDERGARTEN + WITHIN_350M_CHILDCARE +
  WITHIN_350M_BUS + WITHIN_1KM_PRISCH

gwr_predict_at <- function(bw) {
  g <- gwr.basic(f, data = train_data, regression.points = test_data, bw = bw,
                 kernel = "gaussian", adaptive = TRUE, longlat = FALSE)
  B <- as.matrix(st_drop_geometry(st_as_sf(g$SDF))[, 1:15])
  X <- cbind(1, as.matrix(st_drop_geometry(test_data)[, all.vars(f)[-1]]))
  rowSums(B * X)
}

if (!file.exists(p("gwr_pred_bw200.rds"))) {
  t0 <- Sys.time()
  write_rds(tibble(prediction = gwr_predict_at(200)), p("gwr_pred_bw200.rds"))
  cat("GWR bw = 200 prediction:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
pred40 <- read_rds(p("gwr_pred.rds"))$prediction
pred200 <- read_rds(p("gwr_pred_bw200.rds"))$prediction
y <- test_data$resale_price

bandwidth <- tibble(
  bw = c(40, 200),
  RMSE = c(rmse(y, pred40), rmse(y, pred200)),
  MAE = c(mae(y, pred40), mae(y, pred200)),
  MAPE_pct = 100 * c(mape(y, pred40), mape(y, pred200))
)

out <- list(overlap = overlap, twins = twins, bandwidth = bandwidth)
write_rds(out, p("handout_changes.rds"))
print(out)
