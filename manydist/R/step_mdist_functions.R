#' Add manydist dissimilarities to a recipe
#'
#' `step_mdist()` is a [recipes::recipe()] step that replaces selected
#' predictors by a distance-based representation computed with [mdist()]. It is
#' designed for distance-based learning workflows, especially nearest-neighbour
#' prediction and clustering models that operate on dissimilarity matrices.
#'
#' The step can produce either distances from new observations to the training
#' observations, or the within-training pairwise dissimilarity matrix. The
#' former is the usual choice for supervised prediction workflows; the latter
#' is useful for distance-based clustering workflows fitted on the training
#' data.
#'
#' @param recipe A recipe object.
#' @param ... Selector(s) for the predictor columns to be used in [mdist()].
#'   These are passed to [recipes::recipes_eval_select()] during preparation.
#' @param role Role for the new distance columns. The default is `"predictor"`.
#' @param trained Logical for recipes internals. Do not set manually.
#' @param output Character string specifying the type of distance output.
#'   `"distance_to_training"` returns distances from the baked data to the
#'   training observations and is the usual choice for prediction workflows.
#'   `"pairwise"` returns the within-training pairwise dissimilarity matrix and
#'   is intended for training-only distance-based clustering workflows.
#' @param response_used Logical. If `TRUE` (the default) and the recipe has
#'   exactly one outcome, response-aware distance specifications use that
#'   outcome when the step is prepared. Set to `FALSE` to construct the
#'   distance from predictors only. Specifications that are not response-aware
#'   never use the outcome.
#' @param preset Character string specifying the distance preset passed to
#'   [mdist()]. Available values include `"custom"`, `"gower"`,
#'   `"unbiased_dependent"`, `"u_dep"`, `"u_indep"`, `"u_mix"`, `"hl"`,
#'   `"gudmm"`, `"dkss"`, `"mod_gower"`, and `"euclidean"`. Preset parameters
#'   are normally fixed by the selected preset. The exception is `method_num`
#'   for `preset = "euclidean"` when all predictors selected by the step are
#'   numeric.
#' @param method_cat Character string specifying the categorical-variable
#'   dissimilarity passed to [mdist()] when `preset = "custom"`. Common values
#'   include `"matching"` and `"tvd"`. Use [all_dist_method_specs()] to inspect
#'   available methods.
#' @param method_num Character string or `NULL` specifying numerical-variable
#'   preprocessing. Supported values are `"none"` for no preprocessing,
#'   `"std"` for standard-deviation scaling, `"range"` for range scaling,
#'   `"robust"` for inter-quartile-range-based scaling, and `"pc_scores"` for
#'   principal-component score scaling. With `preset = "custom"`, `NULL`
#'   defaults to `"none"`. With `preset = "euclidean"`, the default is `"std"`;
#'   when all predictors selected by the step are numeric, an explicitly
#'   supplied value can override that default.
#' @param commensurable Logical. If `TRUE`, dissimilarities are scaled so that
#'   the average contribution of each variable to the overall distance is equal
#'   to 1, when supported by the selected distance specification.
#' @param ncomp Integer or `NULL`. Number of principal components to retain
#'   when `method_num = "pc_scores"`. If `NULL`, all available components are
#'   used unless `threshold` is supplied and supported by the underlying method.
#' @param threshold Numeric or `NULL`. Optional cumulative variance threshold
#'   used when `method_num = "pc_scores"`.
#' @param columns Names of columns selected at prep time. Used internally by
#'   recipes.
#' @param response_col Name of the outcome selected at prep time when it is
#'   used by the distance specification. Used internally by recipes.
#' @param train_predictors Training predictors stored at prep time. Used
#'   internally by recipes to compute distances from new observations to the
#'   training observations.
#' @param preprocessor Internal fitted manydist preprocessor.
#' @param skip Logical. Standard recipes argument indicating whether the step
#'   should be skipped when baking new data.
#' @param id Character string. Unique step identifier.
#'
#' @details
#' During [recipes::prep()], `step_mdist()` stores the selected training
#' predictors and fits the internal manydist preprocessor. During
#' [recipes::bake()], the selected predictors are removed and replaced by
#' distance columns named `dist_1`, `dist_2`, and so on.
#'
#' When `response_used = TRUE`, the step discovers the outcome from the recipe
#' variable roles. For response-aware presets and custom categorical methods,
#' the outcome from the current analysis set is supplied to [mdist()] during
#' preparation. The fitted response-aware profiles are then reused when new
#' data are baked; assessment and test outcomes are never required or used.
#'
#' With `output = "distance_to_training"`, baking the training data returns the
#' training pairwise distances, while baking new data returns distances from
#' each new observation to each training observation. This rectangular
#' representation is suitable for nearest-neighbour prediction models.
#'
#' With `output = "pairwise"`, the step returns the within-training pairwise
#' dissimilarity matrix. Baking genuinely new data is not supported in this
#' mode, because the output is intended for training-only clustering workflows
#' such as [pam_dist()] or [spectral_dist()].
#'
#' @return An updated recipe with a manydist step.
#'
#' @seealso [mdist()], [nearest_neighbor_dist()], [pam_dist()],
#'   [spectral_dist()], [all_dist_method_specs()]
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
#'   # Distance-to-training representation for prediction workflows
#'   rec <- recipes::recipe(species ~ ., data = penguins_small) |>
#'     step_mdist(
#'       recipes::all_predictors(),
#'       preset = "gower",
#'       output = "distance_to_training"
#'     )
#'
#'   rec_prep <- recipes::prep(rec, training = penguins_small)
#'   baked <- recipes::bake(rec_prep, new_data = penguins_small)
#'
#'   baked |> dplyr::slice_head(n=5)
#'
#'   # Pairwise representation for clustering workflows
#'   rec_pairwise <- recipes::recipe(~ ., data = penguins_small) |>
#'     step_mdist(
#'       recipes::all_predictors(),
#'       preset = "gower",
#'       output = "pairwise"
#'     )
#'
#'   rec_pairwise_prep <- recipes::prep(rec_pairwise, training = penguins_small)
#'   pairwise_dist <- recipes::bake(rec_pairwise_prep, new_data = penguins_small)
#'
#'   pairwise_dist
#' }
#'
#' @export
step_mdist <- function(
    recipe,
    ...,
    role             = "predictor",
    trained          = FALSE,
    output           = "distance_to_training",
    response_used    = TRUE,
    preset           = "custom",
    method_cat       = "tot_var_dist",
    method_num       = NULL,
    commensurable    = FALSE,
    ncomp            = NULL,
    threshold        = NULL,
    columns          = NULL,
    response_col     = NULL,
    train_predictors = NULL,
    preprocessor     = NULL,
    skip             = FALSE,
    id               = recipes::rand_id("mdist")
) {

  commensurable_missing <- missing(commensurable)
  output <- rlang::arg_match0(output, c("distance_to_training", "pairwise"))

  if (!is.logical(response_used) || length(response_used) != 1L ||
      is.na(response_used)) {
    rlang::abort("`response_used` must be a single `TRUE` or `FALSE` value.")
  }

  if (is.null(method_num)) {
    method_num <- if (identical(preset, "custom")) {
      "none"
    } else {
      "std"
    }
  }

  if (commensurable_missing && !identical(preset, "custom")) {
    commensurable <- TRUE
  }

  recipes::add_step(
    recipe,
    step_mdist_new(
      terms            = rlang::enquos(...),
      role             = role,
      trained          = trained,
      output           = output,
      response_used    = response_used,
      preset           = preset,
      method_cat       = method_cat,
      method_num       = method_num,
      commensurable    = commensurable,
      ncomp            = ncomp,
      threshold        = threshold,
      columns          = columns,
      response_col     = response_col,
      train_predictors = train_predictors,
      preprocessor     = preprocessor,
      skip             = skip,
      id               = id
    )
  )
}

step_mdist_new <- function(terms, role, trained,
                           output, response_used,
                           preset, method_cat, method_num,
                           commensurable,
                           ncomp, threshold,
                           columns, response_col,
                           train_predictors, preprocessor,
                           skip, id) {
  recipes::step(
    subclass         = "mdist",
    terms            = terms,
    role             = role,
    trained          = trained,
    output           = output,
    response_used    = response_used,
    preset           = preset,
    method_cat       = method_cat,
    method_num       = method_num,
    commensurable    = commensurable,
    ncomp            = ncomp,
    threshold        = threshold,
    columns          = columns,
    response_col     = response_col,
    train_predictors = train_predictors,
    preprocessor     = preprocessor,
    skip             = skip,
    id               = id
  )
}


.step_mdist_is_response_aware <- function(preset, method_cat) {
  preset <- .normalize_preset(preset)

  if (preset %in% c("unbiased_dependent", "u_dep", "u_mix")) {
    return(TRUE)
  }

  if (!identical(preset, "custom")) {
    return(FALSE)
  }

  method_cat <- ifelse(
    method_cat %in% c("tot_var_dist", "tvd"),
    "tvd",
    method_cat
  )

  any(method_cat %in% response_aware_methods(argument = "method_cat"))
}


#' @export
prep.step_mdist <- function(x, training, info = NULL, ...) {
  col_ids   <- recipes::recipes_eval_select(x$terms, training, info)
  col_names <- names(col_ids)

  if (length(col_names) == 0L) {
    rlang::abort("`step_mdist()` did not select any columns.")
  }

  train_predictors <- training[, col_names, drop = FALSE]
  ncomp <- x$ncomp %||% ncol(train_predictors)
  output <- rlang::arg_match0(x$output, c("distance_to_training", "pairwise"))

  response_col <- NULL
  use_recipe_response <- isTRUE(x$response_used) &&
    .step_mdist_is_response_aware(x$preset, x$method_cat)

  if (use_recipe_response && !is.null(info)) {
    outcome_names <- unique(info$variable[info$role == "outcome"])
    outcome_names <- intersect(outcome_names, colnames(training))

    if (length(outcome_names) > 1L) {
      rlang::abort(
        paste0(
          "Response-aware `step_mdist()` specifications require exactly one ",
          "recipe outcome; found ", length(outcome_names), "."
        )
      )
    }

    if (length(outcome_names) == 1L) {
      response_col <- outcome_names
    }
  }

  if (!is.null(response_col) && response_col %in% col_names) {
    rlang::abort(
      "The recipe outcome selected by `step_mdist()` cannot also be a predictor."
    )
  }

  preprocessor_training <- train_predictors
  if (!is.null(response_col)) {
    preprocessor_training[[response_col]] <- training[[response_col]]
  }

  prep_args <- list(
    x = preprocessor_training,
    method_cat = x$method_cat,
    method_num = x$method_num,
    commensurable = x$commensurable,
    ncomp = ncomp,
    threshold = x$threshold,
    preset = x$preset
  )

  if (!is.null(response_col)) {
    prep_args$response <- response_col
  }

  preprocessor <- rlang::exec(.prep_mdist, !!!prep_args)

  step_mdist_new(
    terms            = x$terms,
    role             = x$role,
    trained          = TRUE,
    output           = output,
    response_used    = x$response_used,
    preset = preprocessor$preset,
    method_cat = preprocessor$method_cat,
    method_num = preprocessor$method_num,
    commensurable = preprocessor$commensurable,
    ncomp            = ncomp,
    threshold        = x$threshold,
    columns          = col_names,
    response_col     = response_col,
    train_predictors = train_predictors,
    preprocessor     = preprocessor,
    skip             = x$skip,
    id               = x$id
  )
}

#' @export
bake.step_mdist <- function(object, new_data, ...) {
  recipes::check_new_data(object$columns, object, new_data)

  pred_mat <- new_data[, object$columns, drop = FALSE]

  same_n     <- nrow(pred_mat) == nrow(object$train_predictors)
  same_p     <- ncol(pred_mat) == ncol(object$train_predictors)
  same_names <- identical(colnames(pred_mat), colnames(object$train_predictors))

  if (same_n && same_p && same_names) {
    is_training <- isTRUE(
      all.equal(
        pred_mat,
        object$train_predictors,
        check.attributes = FALSE
      )
    )
  } else {
    is_training <- FALSE
  }

  if (identical(object$output, "pairwise")) {
    if (!is_training) {
      rlang::abort(
        paste0(
          "`step_mdist(output = 'pairwise')` is intended for within-training ",
          "dissimilarities. Baking new data is not supported directly; ",
          "use the fitted preprocessor internally at prediction time."
        )
      )
    }

    d_mat <- .apply_mdist(object$preprocessor, new_data = NULL)
    d_mat <- as.matrix(d_mat)
    n_train <- nrow(object$train_predictors)
    colnames(d_mat) <- paste0("dist_", seq_len(n_train))

    keep_cols <- setdiff(colnames(new_data), object$columns)

    new_data_out <- cbind(
      new_data[, keep_cols, drop = FALSE],
      as.data.frame(d_mat)
    )

    return(tibble::as_tibble(new_data_out))
  }

  # output == "distance_to_training"
  if (is_training) {
    d_mat <- .apply_mdist(object$preprocessor, new_data = NULL)
  } else {
    d_mat <- .apply_mdist(object$preprocessor, new_data = pred_mat)
  }

  d_mat <- as.matrix(d_mat)
  n_train <- nrow(object$train_predictors)
  colnames(d_mat) <- paste0("dist_", seq_len(n_train))

  keep_cols <- setdiff(colnames(new_data), object$columns)

  new_data_out <- cbind(
    new_data[, keep_cols, drop = FALSE],
    as.data.frame(d_mat)
  )

  tibble::as_tibble(new_data_out)
}

#' @export
print.step_mdist <- function(x, width = max(20, getOption("width") - 30), ...) {

  cat("Step: mdist\n")
  cat("  role:        ", x$role,   "\n", sep = "")
  cat("  trained:     ", x$trained,"\n\n", sep = "")

  cat("  output:      ", x$output, "\n", sep = "")
  cat("  preset:      ", x$preset, "\n", sep = "")

  response_status <- if (!isTRUE(x$response_used)) {
    "<disabled>"
  } else if (!is.null(x$response_col)) {
    x$response_col
  } else if (!isTRUE(x$trained) &&
             .step_mdist_is_response_aware(x$preset, x$method_cat)) {
    "<recipe outcome at prep>"
  } else {
    "<not used>"
  }

  cat("  response:    ", response_status, "\n", sep = "")

  if (x$preset == "custom") {
    cat("  method_cat:    ", x$method_cat,     "\n", sep = "")
    cat("  method_num:    ", x$method_num,     "\n", sep = "")
    cat("  commensurable: ", x$commensurable,  "\n", sep = "")
    cat("  ncomp:         ", x$ncomp,          "\n", sep = "")
    cat("  threshold:     ", x$threshold,      "\n", sep = "")
  } else {
    cat("  (arguments handled internally by preset)\n")
  }

  if (x$trained) {
    cat("\n  columns:        ", paste(x$columns, collapse = ", "), "\n", sep = "")
    cat("  # training rows:", nrow(x$train_predictors), "\n", sep = "")
  }

  invisible(x)
}
