#' Distance-based nearest-neighbour model
#'
#' `nearest_neighbor_dist()` defines a parsnip model specification for
#' nearest-neighbour prediction using precomputed distance representations,
#' such as those produced by [step_mdist()]. It can be used for classification
#' or regression workflows in which predictors have first been transformed into
#' distances to the training observations.
#'
#' This model is intended to be used together with [step_mdist()] in a
#' [recipes::recipe()]. The recipe creates distance columns named `dist_1`,
#' `dist_2`, and so on; `nearest_neighbor_dist()` then applies k-nearest
#' neighbours to that distance representation.
#'
#' @param mode Character string specifying the model mode. Available values are
#'   `"classification"` and `"regression"`.
#' @param neighbors Number of neighbours. This can be an integer or a tunable
#'   parameter, for example `tune::tune()`.
#'
#' @details
#' The model uses a manydist-specific parsnip engine. In the usual workflow,
#' distances are computed by [step_mdist()] with
#' `output = "distance_to_training"`. The resulting distance columns are then
#' passed to `nearest_neighbor_dist()`.
#'
#' Lower-level engine functions such as `fit_knn_dist()` and
#' `predict_knn_dist_*()` are exported for parsnip registration, but users
#' normally do not need to call them directly.
#'
#' @return A parsnip model specification of class `"nearest_neighbor_dist"`.
#'
#' @seealso [step_mdist()], [mdist()]
#'
#' @examples
#' if (requireNamespace("palmerpenguins", quietly = TRUE)) {
#'   data("penguins", package = "palmerpenguins")
#'
#'   penguins_small <- palmerpenguins::penguins |>
#'     dplyr::select(
#'       species, bill_length_mm, bill_depth_mm, flipper_length_mm,
#'       body_mass_g, island, sex
#'     ) |>
#'     tidyr::drop_na()
#'
#'   set.seed(123)
#'   penguin_split <- rsample::initial_split(
#'     penguins_small,
#'     prop = 0.75,
#'     strata = species
#'   )
#'
#'   penguin_train <- rsample::training(penguin_split)
#'   penguin_test  <- rsample::testing(penguin_split)
#'
#'   rec <- recipes::recipe(species ~ ., data = penguin_train) |>
#'     step_mdist(
#'       recipes::all_predictors(),
#'       preset = "gower",
#'       output = "distance_to_training"
#'     )
#'
#'   spec <- nearest_neighbor_dist(
#'     mode = "classification",
#'     neighbors = 5
#'   ) |>
#'     parsnip::set_engine("manydist")
#'
#'   wf <- workflows::workflow() |>
#'     workflows::add_recipe(rec) |>
#'     workflows::add_model(spec)
#'
#'   fit <- workflows::fit(wf, data = penguin_train)
#'
#'   predict(fit, new_data = penguin_test) |>
#'     dplyr::slice_head(n = 5)
#' }
#' @export
nearest_neighbor_dist <- function(mode      = "classification",
                                  neighbors = NULL) {

  if (!mode %in% c("classification", "regression", "unknown")) {
    rlang::abort("`mode` must be 'classification', 'regression', or 'unknown'.")
  }

  args <- list(
    neighbors = rlang::enquo(neighbors)
  )

  parsnip::new_model_spec(
    "nearest_neighbor_dist",
    args     = args,
    eng_args = NULL,
    mode     = mode,
    method   = NULL,
    engine   = "manydist"
  )
}
## =====================================================================
## Parsnip model registration: nearest_neighbor_dist
## =====================================================================
#' @keywords internal
register_nearest_neighbor_dist <- function() {

  existing <- try(parsnip::get_model_env()$models$model,
                  silent = TRUE)

  if (!inherits(existing, "try-error") &&
      "nearest_neighbor_dist" %in% existing) {
    return(invisible(NULL))
  }

  parsnip::set_new_model("nearest_neighbor_dist")

  parsnip::set_model_mode("nearest_neighbor_dist", "classification")
  parsnip::set_model_mode("nearest_neighbor_dist", "regression")

  parsnip::set_model_engine(
    "nearest_neighbor_dist",
    mode = "classification",
    eng  = "manydist"
  )

  parsnip::set_model_engine(
    "nearest_neighbor_dist",
    mode = "regression",
    eng  = "manydist"
  )

  # (no set_dependency here to avoid any side-effects during dev)

  parsnip::set_model_arg(
    model    = "nearest_neighbor_dist",
    eng  = "manydist",
    parsnip  = "neighbors",
    original = "k",          # argument name in fit_knn_dist
    func     = list(pkg = "rlang", fun = "quo"),
    has_submodel = FALSE
  )

  # parsnip::set_model_arg(
  #   model    = "nearest_neighbor_dist",
  #   eng  = "manydist",
  #   parsnip  = "dist_fun",
  #   original = "dist_fun",
  #   func     = list(pkg = "rlang", fun = "quo"),
  #   has_submodel = FALSE
  # )
  #
  # parsnip::set_model_arg(
  #   model    = "nearest_neighbor_dist",
  #   eng  = "manydist",
  #   parsnip  = "dist_args",
  #   original = "dist_args",
  #   func     = list(pkg = "rlang", fun = "quo"),
  #   has_submodel = FALSE
  # )

  parsnip::set_fit(
    model = "nearest_neighbor_dist",
    eng   = "manydist",
    mode  = "classification",
    value = list(
      interface = "data.frame",
      protect   = c("x", "y"),
      func      = c(pkg = "manydist", fun = "fit_knn_dist"),
      defaults  = list()
    )
  )

  parsnip::set_fit(
    model = "nearest_neighbor_dist",
    eng   = "manydist",
    mode  = "regression",
    value = list(
      interface = "data.frame",
      protect   = c("x", "y"),
      func      = c(pkg = "manydist", fun = "fit_knn_dist"),
      defaults  = list()
    )
  )

  parsnip::set_encoding(
    model = "nearest_neighbor_dist",
    eng   = "manydist",
    mode  = "classification",
    options = list(
      predictor_indicators = "none",
      compute_intercept    = FALSE,
      remove_intercept     = FALSE,
      allow_sparse_x       = FALSE
    )
  )

  parsnip::set_encoding(
    model = "nearest_neighbor_dist",
    eng   = "manydist",
    mode  = "regression",
    options = list(
      predictor_indicators = "none",
      compute_intercept    = FALSE,
      remove_intercept     = FALSE,
      allow_sparse_x       = FALSE
    )
  )

  # Debug: inspect the internal registry
  # cat("\n--- internal model registry ---\n")
  # reg <- get("models", envir = parsnip:::get_model_env())
  # str(reg)

  class_info <- list(
    pre  = NULL,
    post = NULL,
    func = c(pkg = "manydist", fun = "predict_knn_dist_class"),
    args = list(
      object   = rlang::expr(object$fit),
      new_data = rlang::expr(new_data),
      type     = "class"
    )
  )

  parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "manydist",
    mode  = "classification",
    type  = "class",
    value = class_info
  )

  prob_info <- list(
    pre  = NULL,
    post = NULL,
    func = c(pkg = "manydist", fun = "predict_knn_dist_prob"),
    args = list(
      object   = rlang::expr(object$fit),
      new_data = rlang::expr(new_data),
      type     = "prob"
    )
  )

  parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "manydist",
    mode  = "classification",
    type  = "prob",
    value = prob_info
  )

  reg_info <- list(
    pre  = NULL,
    post = NULL,
    func = c(pkg = "manydist", fun = "predict_knn_dist_reg"),
    args = list(
      object   = rlang::expr(object$fit),
      new_data = rlang::expr(new_data),
      type     = "numeric"
    )
  )

  parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "manydist",
    mode  = "regression",
    type  = "numeric",
    value = reg_info
  )

}
