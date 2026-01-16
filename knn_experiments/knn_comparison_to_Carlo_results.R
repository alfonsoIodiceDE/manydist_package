## ============================================================
## Compare Carlo's last replication with tidymodels + step_mdist
## Assumes you've already done:
##   load("unbiased_knn_balanced.RData")
## which gives: tr_df, ts_df, y_tr, y_ts, dist_g, dist_oh,
##              dist_udep, dist_uind, acc, k_max, etc.
## ============================================================

rm(list = ls())
load("../knn_experiments/knn_unbiased/unbiased_knn_balanced.RData")

devtools::load_all("../manydist/")

library(tidyverse)
library(tidymodels)

# If registration isn't automatic, uncomment:
# manydist::register_nearest_neighbor_dist()

## ------------------------------------------------------------
## 1. Rebuild train/test data frames with outcome
## ------------------------------------------------------------

train_df_full <- tr_df
train_df_full$y <- y_tr

test_df_full <- ts_df
test_df_full$y <- y_ts

## ------------------------------------------------------------
## 2. Method mapping: Carlo's 4 methods -> mdist settings
## ------------------------------------------------------------

method_map <- tibble::tibble(
  method_code   = c("g",      "oh",        "udep",      "udep",               "uind"),
  dist_obj_name = c("dist_g", "dist_oh",  "dist_udep",  "dist_udep",         "dist_uind"),
  acc_col       = c(1L,       2L,                3L,                  4L, 5L),
  mdist_type    = c("preset", "preset",          "preset",    "custom"  ,      "custom"),
  mdist_preset  = c("gower",  "euclidean_onehot","unbiased_dependent",NA_character_ , NA_character_),
  param_set     = list(
    NULL,  # g
    NULL,  # oh
    NULL,  # udep
    # uind: unbiased independent, Carlo-style
    list(
      preset        = "custom",
      distance_cont = "manhattan",
      distance_cat  = "tot_var_dist",
      scaling_cont="std",
      commensurable = TRUE
      # add scaling_cont etc. here if you need to match old defaults
    ),
    list(
      preset        = "custom",
      distance_cont = "manhattan",
      distance_cat  = "cat_dis",
      commensurable = TRUE
      # add scaling_cont etc. here if you need to match old defaults
    )
  )
)

## Result containers
dist_comp_results <- list()
knn_comp_results  <- list()
acc[]
## ------------------------------------------------------------
## 3. Loop over the four distance methods
## ------------------------------------------------------------

for (ii in seq_len(nrow(method_map))) {


  this_method   <- method_map$method_code[ii]
  this_distname <- method_map$dist_obj_name[ii]
  acc_col       <- method_map$acc_col[ii]
  mdist_type    <- method_map$mdist_type[ii]
  mdist_preset  <- method_map$mdist_preset[ii]
  param_set     <- method_map$param_set[[ii]]  # NULL or list

  cat("\n=============================\n")
  cat("Method:", this_method, "\n")
  cat("=============================\n")

  ## -------------------------------
  ## 3a. Build recipe with step_mdist
  ## -------------------------------
  rec_last <- manydist::make_mdist_recipe(
    df           = train_df_full,
    mdist_type   = mdist_type,
    mdist_preset = mdist_preset,
    param_set    = param_set,
    outcome      = "y"
  )

  prepped_last <- prep(rec_last)
  baked_train  <- bake(prepped_last, new_data = NULL)
  baked_test   <- bake(prepped_last, new_data = test_df_full)

  ## Extract tidymodels distance matrix: rows = test, cols = train
  dist_cols <- grep("^dist_", names(baked_test), value = TRUE)
  dist_test_train_tm <- as.matrix(baked_test[, dist_cols, drop = FALSE])

  cat("tidymodels distance matrix dim (test x train): ",
      paste(dim(dist_test_train_tm), collapse = " x "), "\n")

  ## -------------------------------
  ## 3b. Compare to Carlo's distance
  ## -------------------------------
  dist_carlo <- get(this_distname)  # stored in .RData
  cat("Carlo distance matrix dim (train x test): ",
      paste(dim(dist_carlo), collapse = " x "), "\n")


  # Carlo: train x test; tidymodels: test x train
  # So we compare t(dist_carlo) to dist_test_train_tm
  max_abs_diff_dist <- max(abs(dist_test_train_tm - t(dist_carlo)))
  cat("Max |tidymodels_dist - t(Carlo_dist)| =", max_abs_diff_dist, "\n")

  dist_comp_results[[this_method]] <- list(
    method           = this_method,
    max_abs_diff_dist = max_abs_diff_dist
  )

  ## --------------------------------------------------------
  ## 3c. KNN using Carlo's knn_dist() on tidymodels distances
  ## --------------------------------------------------------
  # Carlo's knn_dist expects: rows = train, cols = test
  dist_for_carlo_tm <- t(dist_test_train_tm)

  k_max_local <- k_max  # from .RData
  acc_k_from_tm_dists <- numeric(k_max_local)

  for (k in 1:k_max_local) {
    preds_k <- knn_dist(
      dist_matrix  = dist_for_carlo_tm,
      train_labels = y_tr,
      k            = k
    )
    acc_k_from_tm_dists[k] <- compute_accuracy(preds_k, y_ts)
  }


  if(acc_col>=4){acc_carlo=acc[,acc_col-1]}else{acc_carlo <- acc[, acc_col]}

  max_abs_diff_acc_carlo_vs_tm_dists <-
    max(abs(acc_carlo - acc_k_from_tm_dists))

  cat("Max |accuracy (Carlo acc) - accuracy (knn_dist on tm dist)| = ",
      max_abs_diff_acc_carlo_vs_tm_dists, "\n")

  ## --------------------------------------------------------
  ## 3d. KNN via tidymodels engine nearest_neighbor_dist
  ## --------------------------------------------------------
  acc_k_tm_engine <- numeric(k_max_local)

  for (k in 1:k_max_local) {
    knn_spec_tm <- nearest_neighbor_dist(
      mode      = "classification",
      neighbors = k
    )

    wf_last <- workflow() |>
      add_recipe(rec_last) |>
      add_model(knn_spec_tm)

    fit_tm  <- fit(wf_last, data = train_df_full)
    pred_tm <- predict(fit_tm, new_data = test_df_full)

    acc_k_tm_engine[k] <- compute_accuracy(pred_tm$.pred_class, y_ts)
  }

  max_abs_diff_acc_carlo_vs_tm_engine <-
    max(abs(acc_carlo - acc_k_tm_engine))

  cat("Max |accuracy (Carlo acc) - accuracy (tidymodels engine)| = ",
      max_abs_diff_acc_carlo_vs_tm_engine, "\n")

  max_abs_diff_acc_dists_vs_engine <-
    max(abs(acc_k_from_tm_dists - acc_k_tm_engine))

  cat("Max |accuracy (knn_dist on tm dist) - accuracy (tm engine)| = ",
      max_abs_diff_acc_dists_vs_engine, "\n")

  knn_comp_results[[this_method]] <- list(
    method                           = this_method,
    acc_carlo                        = acc_carlo,
    acc_from_tm_dists                = acc_k_from_tm_dists,
    acc_from_tm_engine               = acc_k_tm_engine,
    max_abs_diff_acc_carlo_vs_tm_dists  = max_abs_diff_acc_carlo_vs_tm_dists,
    max_abs_diff_acc_carlo_vs_tm_engine = max_abs_diff_acc_carlo_vs_tm_engine,
    max_abs_diff_acc_dists_vs_engine    = max_abs_diff_acc_dists_vs_engine
  )
}

## ------------------------------------------------------------
## 4. Optional: summarize in a tibble
## ------------------------------------------------------------

dist_summary <- map_dfr(dist_comp_results, ~as_tibble(.x))
dist_summary

knn_summary <- map_dfr(knn_comp_results, function(x) {
  tibble(
    method                           = x$method,
    max_abs_diff_acc_carlo_vs_tm_dists  = x$max_abs_diff_acc_carlo_vs_tm_dists,
    max_abs_diff_acc_carlo_vs_tm_engine = x$max_abs_diff_acc_carlo_vs_tm_engine,
    max_abs_diff_acc_dists_vs_engine    = x$max_abs_diff_acc_dists_vs_engine
  )
})
knn_summary
