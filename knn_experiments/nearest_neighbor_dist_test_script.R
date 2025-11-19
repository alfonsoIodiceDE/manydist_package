# AT THE MOMENT THE NN BASED ON A GENERIC DIST FUN SEEMS TO WORK, SO I NEED TO CHECK WHETHER IT STILL WORK
# NEED TO CHECK THE RESPONSE IDENTIFICATION



## =====================================================================
## 0. Packages
## =====================================================================
rm(list = ls())

library(tidymodels)
library(tibble)
library(dplyr)
library(cluster)   # for daisy
library(recipes)   # for gower_recipe + baking
pkgload::load_all(".")
# source("R/knn_dist_functions.R")
# source("R/daisy_gower_dist.R")
# source("R/knn_dist_parsnip_registration.R")


## Dummy distance for dev (train and new_data contain Species here)
dummy_dist_manual <- function(train, new_data, ...) {
  X <- as.matrix(dplyr::select(train, -Species))
  Z <- as.matrix(dplyr::select(new_data, -Species))
  Dfull <- as.matrix(dist(rbind(Z, X)))  # (n_test + n_train) × (n_test + n_train)
  Dfull[seq_len(nrow(Z)), (nrow(Z) + 1):(nrow(Z) + nrow(X))]  # n_test × n_train
}

dummy_dist_parsnip <- function(train, new_data, ...) {
  X <- as.matrix(train)
  Z <- as.matrix(new_data)
  Dfull <- as.matrix(dist(rbind(Z, X)))
  Dfull[seq_len(nrow(Z)), (nrow(Z) + 1):(nrow(Z) + nrow(X))]
}
## Quick sanity check of base code -------------------------------------

dat <- as_tibble(iris) %>% mutate(Species = factor(Species))
set.seed(1)
sp  <- initial_split(dat, prop = 0.7, strata = Species)
data_tr  <- training(sp)
data_ts  <- testing(sp)



fit_daisy <- fit_knn_dist(
  x         = data_tr,
  y         = data_tr$Species,
  k         = 5,
  dist_fun  = daisy_gower_dist,
  dist_args = list()
)

pred_daisy <- predict_knn_dist_class(fit_daisy, data_ts, type = "class")

acc_daisy <- yardstick::accuracy_vec(
  truth   = data_ts$Species,
  estimate = pred_daisy
)
acc_daisy

mdist_gower_dist <- function(train, new_data, ...) {
  X <- dplyr::select(train, -Species)
  Z <- dplyr::select(new_data, -Species)

  # your mdist must have signature mdist(x, new_data = NULL, preset = ...)
  md <- mdist(x = X, new_data = Z, preset = "gower")

  # return the MDist object; predict_knn_dist_class knows how to handle it
  md
}

fit_md <- fit_knn_dist(
  x         = data_tr,
  y         = data_tr$Species,
  k         = 5,
  dist_fun  = mdist_gower_dist,
  dist_args = list()
)



pred_md <- predict_knn_dist_class(fit_md, data_ts, type = "class")

acc_md <- yardstick::accuracy_vec(
  truth   = data_ts$Species,
  estimate = pred_md
)
acc_md

fit0 <- fit_knn_dist(x = tr, y = tr$Species, k = 5,
                     dist_fun = dummy_dist_manual, dist_args = list())

pred_class0 <- predict_knn_dist_class(fit0, te, type = "class")
acc0 <- yardstick::accuracy_vec(truth = te$Species, estimate = pred_class0)
cat("Base dummy kNN accuracy:", acc0, "\n\n")




# iris data
dat <- as_tibble(iris) %>% mutate(Species = factor(Species))
set.seed(1)
sp  <- initial_split(dat, prop = 0.7, strata = Species)
tr  <- training(sp)
te  <- testing(sp)


spec <- nearest_neighbor_dist(
  mode      = "classification",
  neighbors = 5,
  dist_fun  = dummy_dist_parsnip,
  dist_args = list()
) |>
  set_engine("precomputed")

fit_parsnip <- fit(
  spec,
  Species ~ .,
  data = tr
)

str(fit_parsnip$fit)

pred_class <- predict(fit_parsnip, new_data = te, type = "class")

show_model_info("nearest_neighbor_dist")



yardstick::accuracy(
  bind_cols(te |> select(Species), pred_class),
  truth   = Species,
  estimate = .pred_class
)

pkgload::load_all(".")

mdist_parsnip <- function(train, new_data, ...) {
  # train and new_data: predictors only
  # call your local mdist, making sure it uses
  #   mdist(x = train, new_data = new_data, preset = ..., ...)
  md <- mdist(x = train,new_data = new_data,...)
  # md is an MDist, return that (predict_* already knows how to handle MDist)
  md
}

spec_m <- nearest_neighbor_dist(
  mode      = "classification",
  neighbors = 5,
  dist_fun  = mdist_parsnip,
  dist_args = list(preset = "gower")  # or whatever
) |>
  set_engine("precomputed")

fit_m <- fit(spec_m, Species ~ ., data = tr)

pred_m <- predict(fit_m, new_data = te, type = "class")
table(pred_m$.pred_class)

accuracy(bind_cols(te |> select(Species), pred_m),
         truth = Species, estimate = .pred_class)
