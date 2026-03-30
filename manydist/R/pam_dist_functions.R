# R/pam_dist_functions.R

# internal helpers ---------------------------------------------------------

.get_step_mdist <- function(prepped_recipe) {
  idx <- purrr::detect_index(prepped_recipe$steps, ~ inherits(.x, "step_mdist"))
  if (idx == 0) {
    stop("No `step_mdist()` found in the recipe.", call. = FALSE)
  }
  prepped_recipe$steps[[idx]]
}

.extract_dist_block <- function(x, prefix = "^dist_") {
  x <- tibble::as_tibble(x)
  dist_cols <- grep(prefix, names(x), value = TRUE)

  if (length(dist_cols) == 0L) {
    stop("No distance columns found in baked data.", call. = FALSE)
  }

  as.matrix(x[, dist_cols, drop = FALSE])
}

.fit_pam_from_recipe <- function(spec, recipe, data) {
  if (missing(recipe) || missing(data)) {
    stop("`recipe` and `data` must be provided.", call. = FALSE)
  }

  prepped_recipe <- recipes::prep(recipe, training = data)
  mdist_step <- .get_step_mdist(prepped_recipe)

  if (!identical(mdist_step$output, "pairwise")) {
    stop(
      "`pam_dist` requires `step_mdist(..., output = 'pairwise')` in the recipe.",
      call. = FALSE
    )
  }

  baked_train <- recipes::bake(prepped_recipe, new_data = data)
  D_train <- .extract_dist_block(baked_train)

  if (nrow(D_train) != ncol(D_train)) {
    stop("The baked training dissimilarity matrix is not square.", call. = FALSE)
  }

  k_val <- rlang::eval_tidy(spec$num_clusters)

  if (!is.numeric(k_val) || length(k_val) != 1L || is.na(k_val) || k_val < 2) {
    stop("`num_clusters` must evaluate to a single integer >= 2.", call. = FALSE)
  }

  k_val <- as.integer(k_val)

  pam_fit <- cluster::pam(
    x = stats::as.dist(D_train),
    k = k_val,
    diss = TRUE
  )

  structure(
    list(
      spec           = spec,
      fit            = pam_fit,
      recipe         = recipe,
      prepped_recipe = prepped_recipe,
      mdist_step     = mdist_step,
      train_matrix   = D_train,
      medoid_idx     = pam_fit$id.med,
      clusters       = pam_fit$clustering,
      num_clusters   = k_val
    ),
    class = "pam_dist_fit"
  )
}


# specification ------------------------------------------------------------

#' PAM clustering specification based on manydist dissimilarities
#'
#' @param num_clusters Number of clusters.
#'
#' @return A `pam_dist_spec` object.
#' @export
pam_dist <- function(num_clusters = NULL) {
  structure(
    list(
      num_clusters = rlang::enquo(num_clusters),
      engine = "manydist"
    ),
    class = "pam_dist_spec"
  )
}

#' @exportS3Method base::print
print.pam_dist_spec <- function(x, ...) {
  k_val <- tryCatch(rlang::eval_tidy(x$num_clusters), error = function(e) NULL)

  cat("PAM distance clustering specification\n")
  cat("  num_clusters :", if (is.null(k_val)) "<tune/expr>" else k_val, "\n")
  cat("  engine       :", x$engine, "\n")
  invisible(x)
}

#' @exportS3Method parsnip::set_engine
set_engine.pam_dist_spec <- function(object, engine, ...) {
  object$engine <- engine
  object
}


# fit ----------------------------------------------------------------------

#' @exportS3Method generics::fit
fit.pam_dist_spec <- function(object, x, y = NULL, recipe = NULL, data = NULL, ...) {
  engine <- object$engine %||% "manydist"

  if (!identical(engine, "manydist")) {
    stop("Only `engine = 'manydist'` is currently supported.", call. = FALSE)
  }

  if (!is.null(recipe) && !is.null(data)) {
    return(.fit_pam_from_recipe(
      spec = object,
      recipe = recipe,
      data = data
    ))
  }

  if (inherits(x, "recipe")) {
    if (is.null(data)) {
      stop("When `x` is a recipe, `data` must be supplied.", call. = FALSE)
    }
    return(.fit_pam_from_recipe(
      spec = object,
      recipe = x,
      data = data
    ))
  }

  stop(
    paste0(
      "`generics::fit()` for `pam_dist` expects either:\n",
      "  - `fit(spec, recipe = <recipe>, data = <data>)`, or\n",
      "  - `fit(spec, x = <recipe>, data = <data>)`."
    ),
    call. = FALSE
  )
}


# fitted object methods ----------------------------------------------------

#' @exportS3Method base::print
print.pam_dist_fit <- function(x, ...) {
  cat("pam_dist fit\n")
  cat("  num_clusters :", x$num_clusters, "\n")
  cat("  n_obs        :", nrow(x$train_matrix), "\n")
  cat("  medoid idx   :", paste(x$medoid_idx, collapse = ", "), "\n")
  invisible(x)
}

#' @exportS3Method stats::predict
predict.pam_dist_fit <- function(object, new_data = NULL, type = c("cluster", "dist"), ...) {
  type <- match.arg(type)

  if (is.null(new_data)) {
    if (identical(type, "cluster")) {
      return(tibble::tibble(.pred_cluster = factor(object$clusters)))
    }

    D_to_medoids <- object$train_matrix[, object$medoid_idx, drop = FALSE]
    colnames(D_to_medoids) <- paste0("medoid_", seq_len(ncol(D_to_medoids)))
    return(tibble::as_tibble(D_to_medoids))
  }

  step_obj <- object$mdist_step
  pred_mat <- tibble::as_tibble(new_data)[, step_obj$columns, drop = FALSE]

  D_test_train <- .apply_mdist(step_obj$preprocessor, new_data = pred_mat)
  D_test_train <- as.matrix(D_test_train)

  if (ncol(D_test_train) < max(object$medoid_idx)) {
    stop(
      "The computed test-to-training dissimilarity matrix does not contain ",
      "enough columns to identify the fitted medoids.",
      call. = FALSE
    )
  }

  D_test_medoids <- D_test_train[, object$medoid_idx, drop = FALSE]
  colnames(D_test_medoids) <- paste0("medoid_", seq_len(ncol(D_test_medoids)))

  if (identical(type, "dist")) {
    return(tibble::as_tibble(D_test_medoids))
  }

  assigned <- apply(D_test_medoids, 1, which.min)

  tibble::tibble(
    .pred_cluster = factor(assigned, levels = seq_len(object$num_clusters))
  )
}
