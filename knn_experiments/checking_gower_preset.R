library(recipes)
library(cluster)   # for daisy()
library(dplyr)

## ------------------------------------------------------------
## 0. Example continuous-only data
## ------------------------------------------------------------
# You can replace this with your own X
X <- iris %>%
  select(where(is.numeric))   # Sepal.Length, Sepal.Width, ...

## ------------------------------------------------------------
## 1. Apply your gower_recipe and bake the data
## ------------------------------------------------------------
# NOTE: adapt this block to your actual gower_recipe signature.
# I'm assuming gower_recipe() returns a *recipe* object.

# e.g. if gower_recipe expects a data frame of predictors:
rec_gower <- gower_recipe(X)      # <-- adjust arguments if needed
rec_gower_prep <- prep(rec_gower, training = X)
X_baked <- bake(rec_gower_prep, new_data = NULL)

# Sanity check
dim(X)
dim(X_baked)
head(X_baked)

## ------------------------------------------------------------
## 2. Manhattan distance on baked data
## ------------------------------------------------------------
D_manh <- as.matrix(dist(as.matrix(X_baked), method = "manhattan"))

## ------------------------------------------------------------
## 3. Gower distance via daisy() on original data
## ------------------------------------------------------------
D_daisy <- as.matrix(cluster::daisy(X, metric = "gower"))

## ------------------------------------------------------------
## 4. Compare the two distance matrices
## ------------------------------------------------------------
# They should match up to small numerical tolerance
(D_manh/D_daisy)[1:4,1:4]


################################################################################################################
################################################################################################################
################################################################################################################

set.seed(1)
n <- nrow(X)
idx <- sample(n, size = floor(0.7 * n))
X_tr <- X[idx, , drop = FALSE]
X_te <- X[-idx, , drop = FALSE]

## Train: prep on X_tr
rec_gower_tr <- gower_recipe(X_tr)  # adjust signature
rec_gower_tr_prep <- prep(rec_gower_tr, training = X_tr)

X_tr_baked <- bake(rec_gower_tr_prep, new_data = X_tr)
X_te_baked <- bake(rec_gower_tr_prep, new_data = X_te)

# Manhattan test–train:
D_manh_tt <- as.matrix(dist(rbind(X_te_baked, X_tr_baked), method = "manhattan"))
n_te <- nrow(X_te_baked)
n_tr <- nrow(X_tr_baked)
D_manh_rect <- D_manh_tt[1:n_te, (n_te + 1):(n_te + n_tr), drop = FALSE]

# Gower test–train via daisy:
full_orig <- rbind(X_te, X_tr)
D_full_gow <- as.matrix(cluster::daisy(full_orig, metric = "gower"))
D_gow_rect <- D_full_gow[1:n_te, (n_te + 1):(n_te + n_tr), drop = FALSE]

D_manh_rect[1:5,1:5]/D_gow_rect[1:5,1:5]
# Compare
diff_vec_tt <- as.numeric(D_manh_rect - D_gow_rect)
summary(diff_vec_tt)
max(abs(diff_vec_tt), na.rm = TRUE)

################################################################################################################
################################################################################################################
################################################################################################################


# summary(diff_vec)
# max_abs_diff <- max(abs(diff_vec), na.rm = TRUE)
# max_abs_diff
#
#
# X_tr <- tr |> select(-Species)
# X_te <- te |> select(-Species)
#
# md_test <- mdist_parsnip(train = X_tr, new_data = X_te, preset = "gower")
#
# D <- if (inherits(md_test, "MDist")) {
#   as.matrix(md_test$distance)
# } else {
#   as.matrix(md_test)
# }
#
# dim(D)
# # should be 45 x 105
#
# # Are all rows identical?
# all_rows_identical <- all(apply(D, 1, function(row) all(row == D[1, ])))
# all_rows_identical
#
# X_tr_small <- X_tr[1:5, ]
# X_te_small <- X_tr[6:7, ]   # use two rows that are clearly different
#
# md_small <- mdist(x = X_tr_small, new_data = X_te_small, preset = "gower")
# D_small  <- as.matrix(md_small$distance)
#
# D_small
#
#
#
# daisy_gower_parsnip <- function(train, new_data, ...) {
#   require(cluster)
#
#   # full combined dataset
#   full <- rbind(new_data, train)
#
#   # gower dissimilarities
#   Dfull <- as.matrix(cluster::daisy(full, metric = "gower", ...))
#
#   n_test  <- nrow(new_data)
#   n_train <- nrow(train)
#
#   # extract test × train block
#   Dfull[1:n_test, (n_test + 1):(n_test + n_train)]
# }
#
# spec_daisy <- nearest_neighbor_dist(
#   mode      = "classification",
#   neighbors = 5,
#   dist_fun  = daisy_gower_parsnip,
#   dist_args = list()
# ) %>%
#   set_engine("precomputed")
#
# fit_daisy <- fit(spec_daisy, Species ~ ., data = tr)
#
# pred_daisy <- predict(fit_daisy, new_data = te, type = "class")
#
# yardstick::accuracy(
#   tibble(Species = te$Species, .pred_class = pred_daisy$.pred_class),
#   truth = Species,
#   estimate = .pred_class
# )
#
#
# source("R/mdist.R")   # loads your modified mdist() into .GlobalEnv
# pkgload::load_all(".", export_all = TRUE)   # loads package, but won't overwrite the global mdist()
#
# X_tr <- tr |> select(-Species)
# X_te <- te |> select(-Species)
#
# # mdist
# md_m <- mdist(x = X_tr, new_data = X_te, preset = "gower")
# D_m  <- as.matrix(md_m$distance)
# D_m[1:4,1:4]
#
# # daisy (reference)
# all  <- rbind(X_te, X_tr)
# g    <- as.matrix(cluster::daisy(all, metric = "gower"))
# n_te <- nrow(X_te)
# n_tr <- nrow(X_tr)
# D_ref <- g[1:n_te, (n_te + 1):(n_te + n_tr), drop = FALSE]
#
# summary(as.numeric(D_m - D_ref))


mdist_parsnip_core <- function(train, new_data, ...) {
  # predictors only
  md <- mdist(x = train, new_data = new_data, ...)
  md
}

fit_m_core <- fit_knn_dist(
  x         = X_tr,
  y         = tr$Species,
  k         = 5,
  dist_fun  = mdist_parsnip_core,
  dist_args = list(preset = "gower")
)

pred_m_core <- predict_knn_dist_class(fit_m_core, X_te, type = "class")
acc_m_core  <- accuracy_vec(truth = te$Species, estimate = pred_m_core)
acc_m_core
