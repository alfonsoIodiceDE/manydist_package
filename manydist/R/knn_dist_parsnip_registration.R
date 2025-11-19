#' Nearest neighbor with precomputed distances
#'
#' @param mode "classification" or "regression".
#' @param neighbors Number of neighbors (k).
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
    eng_args = NULL,        # engine args come from set_engine()
    mode     = mode,
    method   = NULL,
    engine   = "precomputed"
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
    eng  = "precomputed"
  )

  parsnip::set_model_engine(
    "nearest_neighbor_dist",
    mode = "regression",
    eng  = "precomputed"
  )

  # (no set_dependency here to avoid any side-effects during dev)

  parsnip::set_model_arg(
    model    = "nearest_neighbor_dist",
    eng      = "precomputed",
    parsnip  = "neighbors",
    original = "k",          # argument name in fit_knn_dist
    func     = list(pkg = "rlang", fun = "quo"),
    has_submodel = FALSE
  )

  # parsnip::set_model_arg(
  #   model    = "nearest_neighbor_dist",
  #   eng      = "precomputed",
  #   parsnip  = "dist_fun",
  #   original = "dist_fun",
  #   func     = list(pkg = "rlang", fun = "quo"),
  #   has_submodel = FALSE
  # )
  #
  # parsnip::set_model_arg(
  #   model    = "nearest_neighbor_dist",
  #   eng      = "precomputed",
  #   parsnip  = "dist_args",
  #   original = "dist_args",
  #   func     = list(pkg = "rlang", fun = "quo"),
  #   has_submodel = FALSE
  # )

  parsnip::set_fit(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
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
    eng   = "precomputed",
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
    eng   = "precomputed",
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
    eng   = "precomputed",
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
    eng   = "precomputed",
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
    eng   = "precomputed",
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
    eng   = "precomputed",
    mode  = "regression",
    type  = "numeric",
    value = reg_info
  )

}
