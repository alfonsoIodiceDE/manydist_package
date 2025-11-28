#' Distance-based representation via manydist
#'
#' @param recipe A recipe object.
#' @param ... Selector(s) for the predictor columns to be used in mdist.
#' @param role Role for the new distance columns (usually "predictor").
#' @param trained Logical for recipes internals (do not set manually).
#' @param preset,distance_cont,distance_cat,commensurable,scaling_cont,
#'   ncomp,threshold Arguments forwarded to `mdist()`.
#' @param columns Names of columns selected at prep time (recipes internals).
#' @param train_predictors Training predictors stored at prep time (internals).
#' @param skip,id Standard recipes arguments.
#'
#' @export
step_mdist <- function(
    recipe,
    ...,
    role            = "predictor",
    trained         = FALSE,
    preset          = "custom",
    distance_cont   = "manhattan",
    distance_cat    = "tot_var_dist",
    commensurable   = FALSE,
    scaling_cont    = "none",
    ncomp           = NULL,
    threshold       = NULL,
    columns         = NULL,
    train_predictors = NULL,
    skip            = FALSE,
    id              = recipes::rand_id("mdist")
) {
  recipes::add_step(
    recipe,
    step_mdist_new(
      terms           = rlang::enquos(...),
      role            = role,
      trained         = trained,
      preset          = preset,
      distance_cont   = distance_cont,
      distance_cat    = distance_cat,
      commensurable   = commensurable,
      scaling_cont    = scaling_cont,
      ncomp           = ncomp,
      threshold       = threshold,
      columns         = columns,
      train_predictors = train_predictors,
      skip            = skip,
      id              = id
    )
  )
}

step_mdist_new <- function(terms, role, trained,
                           preset, distance_cont, distance_cat,
                           commensurable, scaling_cont,
                           ncomp, threshold,
                           columns, train_predictors,
                           skip, id) {
  recipes::step(
    subclass        = "mdist",
    terms           = terms,
    role            = role,
    trained         = trained,
    preset          = preset,
    distance_cont   = distance_cont,
    distance_cat    = distance_cat,
    commensurable   = commensurable,
    scaling_cont    = scaling_cont,
    ncomp           = ncomp,
    threshold       = threshold,
    columns         = columns,
    train_predictors = train_predictors,
    skip            = skip,
    id              = id
  )
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

  step_mdist_new(
    terms            = x$terms,
    role             = x$role,
    trained          = TRUE,
    preset           = x$preset,
    distance_cont    = x$distance_cont,
    distance_cat     = x$distance_cat,
    commensurable    = x$commensurable,
    scaling_cont     = x$scaling_cont,
    ncomp            = ncomp,
    threshold        = x$threshold,
    columns          = col_names,
    train_predictors = train_predictors,
    skip             = x$skip,
    id               = x$id
  )
}

#' @export
bake.step_mdist <- function(object, new_data, ...) {

  recipes::check_new_data(object$columns, object, new_data)

  pred_mat <- new_data[, object$columns, drop = FALSE]

  ## ---- Detect if we're baking the training set ----
  same_n     <- nrow(pred_mat) == nrow(object$train_predictors)
  same_p     <- ncol(pred_mat) == ncol(object$train_predictors)
  same_names <- identical(colnames(pred_mat), colnames(object$train_predictors))

  if (same_n && same_p && same_names) {
    # only now we try equality; no warnings on different sizes
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
  ## -------------------------------------------------

  if (object$preset == "custom") {

    if (is_training) {
      # training: within-sample branch in mdist (validate_x is NULL)
      d_obj <- mdist(
        x             = object$train_predictors,
        distance_cont = object$distance_cont,
        distance_cat  = object$distance_cat,
        commensurable = object$commensurable,
        scaling_cont  = object$scaling_cont,
        ncomp         = object$ncomp,
        threshold     = object$threshold,
        preset        = "custom"
      )
    } else {
      # test / new data: train–test distances
      d_obj <- mdist(
        x             = object$train_predictors,
        new_data      = pred_mat,
        distance_cont = object$distance_cont,
        distance_cat  = object$distance_cat,
        commensurable = object$commensurable,
        scaling_cont  = object$scaling_cont,
        ncomp         = object$ncomp,
        threshold     = object$threshold,
        preset        = "custom"
      )
    }

  } else { # preset != "custom"

    if (is_training) {
      d_obj <- mdist(
        x      = object$train_predictors,
        preset = object$preset
      )
    } else {
      d_obj <- mdist(
        x        = object$train_predictors,
        new_data = pred_mat,
        preset   = object$preset
      )
    }

  }

  d_mat <- as.matrix(d_obj$distance)
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

  cat("  preset:      ", x$preset, "\n", sep = "")

  if (x$preset == "custom") {
    # Only show arguments that matter for custom mode
    cat("  distance_cont: ", x$distance_cont,  "\n", sep = "")
    cat("  distance_cat:  ", x$distance_cat,   "\n", sep = "")
    cat("  commensurable: ", x$commensurable,  "\n", sep = "")
    cat("  scaling_cont:  ", x$scaling_cont,   "\n", sep = "")
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
