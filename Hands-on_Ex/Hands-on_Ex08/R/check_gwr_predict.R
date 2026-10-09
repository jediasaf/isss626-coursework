# ---------------------------------------------------------------------------
# check_gwr_predict.R - is the direct GWR prediction the same as gwr.predict()?
#
# On a subset small enough for gwr.predict() to run (3,000 training flats,
# 400 test flats), compare its predictions with local coefficients from
# gwr.basic(regression.points = ...) applied to each test flat's attributes.
# Output: data/model/gwr_predict_check.rds
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(sf); library(GWmodel); library(tidyverse)})
tr <- read_rds("data/model/train_data.rds"); te <- read_rds("data/model/test_data.rds")
f <- formula(read_rds("data/model/price_mlr.rds"))
set.seed(1); a <- tr[sample(nrow(tr), 3000), ]; b <- te[sample(nrow(te), 400), ]

via_package <- gwr.predict(f, data = a, predictdata = b, bw = 40, kernel = "gaussian",
                           adaptive = TRUE, longlat = FALSE)$SDF$prediction
g <- gwr.basic(f, data = a, regression.points = b, bw = 40, kernel = "gaussian",
               adaptive = TRUE, longlat = FALSE)
B <- as.matrix(st_drop_geometry(st_as_sf(g$SDF))[, 1:15])
X <- cbind(1, as.matrix(st_drop_geometry(b)[, all.vars(f)[-1]]))
direct <- rowSums(B * X)

res <- tibble(n_train = nrow(a), n_test = nrow(b),
              max_abs_diff = max(abs(via_package - direct)),
              correlation = cor(via_package, direct))
write_rds(res, "data/model/gwr_predict_check.rds"); print(res)
