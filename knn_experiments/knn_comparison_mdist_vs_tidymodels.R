# Make sure manydist is loaded and model is registered
rm(list=ls())
# load("../knn_experiments/knn_unbiased/unbiased_knn_balanced.RData")
# ls()
devtools::load_all("../manydist/")
library(tidyverse)
library(tidymodels)
library(palmerpenguins)

penguins_knn_mdist_method_check <- function(
    method,
    k         = 5,
    seed      = 123,
    idx_block = 1:5   # rows/cols used for distance sanity check
) {

  # ----------------------------------------------------------
  # 0) Look up mdist settings for this method
  # ----------------------------------------------------------
  row <- mdist_method_lookup |>
    dplyr::filter(method == !!method) |>
    dplyr::slice(1)

  if (nrow(row) == 0) {
    stop("Method '", method, "' not found in mdist_method_lookup.")
  }

  mdist_type   <- row$mdist_type
  mdist_preset <- row$mdist_preset
  param_set    <- row$param_set[[1]]  # NULL or a list

  # ----------------------------------------------------------
  # 1) Data and split (penguins)
  # ----------------------------------------------------------
  set.seed(seed)

  penguins_df <- penguins |>
    dplyr::select(
      species,
      bill_length_mm,
      bill_depth_mm,
      flipper_length_mm,
      body_mass_g,
      island,
      sex
    ) |>
    tidyr::drop_na()

  split <- initial_split(penguins_df, prop = 0.8, strata = species)
  train_data <- training(split)
  test_data  <- testing(split)

  x_train <- dplyr::select(train_data, -species)
  y_train <- train_data$species
  x_test  <- dplyr::select(test_data,  -species)

  # ----------------------------------------------------------
  # 2) Direct mdist() + fit_knn_dist()
  # ----------------------------------------------------------
  # Build mdist arg list

  # Start from param_set if present, else empty list
  if (is.null(param_set)) {
    param_set_md <- list()
  } else {
    param_set_md <- param_set
  }

  # If type is "preset", we set the preset name
  if (mdist_type == "preset") {
    param_set_md$preset <- mdist_preset
  }

  # Add x and new_data for this run
  param_set_md$x        <- x_train
  param_set_md$new_data <- x_test

  # mdist call
  mdist_out <- do.call(mdist, args = param_set_md)
  dist_mat_test <- as.matrix(mdist_out$distance)  # n_test x n_train

  # dist_args for fit_knn_dist: remove x / new_data
  param_set_knn <- param_set_md
  param_set_knn$x        <- NULL
  param_set_knn$new_data <- NULL

  knn_direct_fit <- fit_knn_dist(
    x         = x_train,
    y         = y_train,
    k         = k,
    dist_fun  = mdist,
    dist_args = param_set_knn
  )

  pred_direct <- predict_knn_dist_class(knn_direct_fit, new_data = x_test)
  pred_direct <- factor(pred_direct, levels = levels(y_train))

  # ----------------------------------------------------------
  # 3) Tidymodels: make_mdist_recipe + nearest_neighbor_dist()
  # ----------------------------------------------------------
  rec_penguins <- manydist::make_mdist_recipe(
    df           = train_data,
    mdist_type   = mdist_type,
    mdist_preset = mdist_preset,
    param_set    = param_set,
    outcome      = "species"
  )

  prepped_penguins <- prep(rec_penguins)
  baked_train      <- bake(prepped_penguins, new_data = NULL)
  baked_test       <- bake(prepped_penguins, new_data = test_data)

  # ----------------------------------------------------------
  # 4) Distance sanity check (step_mdist vs direct mdist)
  # ----------------------------------------------------------
  # assume columns dist_1, dist_2, ..., dist_{n_train}
  max_idx_possible <- min(nrow(baked_test), ncol(dist_mat_test))
  idx <- idx_block[idx_block <= max_idx_possible]

  dist_from_step <- baked_test[idx, paste0("dist_", idx), drop = FALSE] |>
    as.matrix()
  dist_from_mdist <- dist_mat_test[idx, idx, drop = FALSE]

  ratio <- dist_from_step / dist_from_mdist
  max_rel_err <- max(abs(ratio - 1), na.rm = TRUE)

  # ----------------------------------------------------------
  # 5) KNN via tidymodels: fit() + predict()
  # ----------------------------------------------------------
  knn_spec_tm <- nearest_neighbor_dist(
    mode      = "classification",
    neighbors = k
  )

  wf_penguins <- workflow() |>
    add_recipe(rec_penguins) |>
    add_model(knn_spec_tm)

  fit_tm  <- fit(wf_penguins, data = train_data)
  pred_tm <- predict(fit_tm, new_data = test_data)

  # ----------------------------------------------------------
  # 6) KNN via tidymodels: last_fit()
  # ----------------------------------------------------------
  last_fit_obj <- last_fit(wf_penguins, split)

  # predictions on the test portion of 'split'
  pred_last_tbl <- collect_predictions(last_fit_obj)
  # assume same row order as `test_data`
  pred_last <- pred_last_tbl$.pred_class

  # ----------------------------------------------------------
  # 7) Compare predictions
  # ----------------------------------------------------------
  comparison <- tibble(
    direct        = pred_direct,
    tidymodels_fit  = pred_tm$.pred_class,
    tidymodels_last = pred_last
  ) |>
    mutate(
      equal_fit_direct  = tidymodels_fit  == direct,
      equal_last_direct = tidymodels_last == direct,
      equal_fit_last    = tidymodels_fit  == tidymodels_last
    )

  summary_tbl <- comparison |>
    summarise(
      prop_equal_fit_direct  = mean(equal_fit_direct),
      n_diff_fit_direct      = sum(!equal_fit_direct),
      prop_equal_last_direct = mean(equal_last_direct),
      n_diff_last_direct     = sum(!equal_last_direct),
      prop_equal_fit_last    = mean(equal_fit_last),
      n_diff_fit_last        = sum(!equal_fit_last),
      n_test                 = dplyr::n(),
      .groups                = "drop_last"
    )

  tibble(
    method               = method,
    mdist_type           = mdist_type,
    mdist_preset         = mdist_preset,
    k                    = k,
    seed                 = seed,
    n_train              = nrow(train_data),
    n_test               = summary_tbl$n_test,
    # direct vs fit()
    prop_equal_fit_direct  = summary_tbl$prop_equal_fit_direct,
    n_diff_fit_direct      = summary_tbl$n_diff_fit_direct,
    # direct vs last_fit()
    prop_equal_last_direct = summary_tbl$prop_equal_last_direct,
    n_diff_last_direct     = summary_tbl$n_diff_last_direct,
    # fit() vs last_fit()
    prop_equal_fit_last    = summary_tbl$prop_equal_fit_last,
    n_diff_fit_last        = summary_tbl$n_diff_fit_last,
    # distance check
    max_rel_err_dist       = max_rel_err,
    param_set              = list(param_set)
  )
}

# Quick checks:
all_methods = c( "naive","gower","hl","hl_add","u_dep","u_ind","u_dep_manh","u_dep_eucl","euclidean_onehot")
  # "gudmm",
  # "dkss",
  # "mod_gower",

tibble(preset_method = all_methods) |> mutate(results=map(.x=preset_method,~penguins_knn_mdist_method_check(.x))) |> unnest(results)

