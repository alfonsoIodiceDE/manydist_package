# internal helpers ---------------------------------------------------------

# .get_step_mdist <- function(prepped_recipe) {
#   idx <- purrr::detect_index(prepped_recipe$steps, ~ inherits(.x, "step_mdist"))
#   if (idx == 0) {
#     stop("No `step_mdist()` found in the recipe.", call. = FALSE)
#   }
#   prepped_recipe$steps[[idx]]
# }
#
# .extract_dist_block <- function(x, prefix = "^dist_") {
#   x <- tibble::as_tibble(x)
#   dist_cols <- grep(prefix, names(x), value = TRUE)
#
#   if (length(dist_cols) == 0L) {
#     stop("No distance columns found in baked data.", call. = FALSE)
#   }
#
#   as.matrix(x[, dist_cols, drop = FALSE])
# }

.fit_pam_dist <- function(spec, recipe, data) {
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

  pam_fit <- cluster::pam(
    x = stats::as.dist(D_train),
    k = spec$num_clusters,
    diss = TRUE
  )

  structure(
    list(
      spec           = spec,
      fit            = pam_fit,
      recipe         = recipe,
      prepped_recipe = prepped_recipe,
      mdist_step     = mdist_step,
      train_data     = tibble::as_tibble(data),
      train_matrix   = D_train,
      medoid_idx     = pam_fit$id.med,
      clusters       = pam_fit$clustering,
      num_clusters   = spec$num_clusters
    ),
    class = "pam_dist_fit"
  )
}


# specification ------------------------------------------------------------

#' PAM clustering specification for manydist dissimilarities
#'
#' @param num_clusters Number of clusters.
#'
#' @return A `pam_dist` model specification.
#' @export
pam_dist <- function(num_clusters) {
  if (missing(num_clusters)) {
    stop("`num_clusters` must be provided.", call. = FALSE)
  }

  if (!is.numeric(num_clusters) || length(num_clusters) != 1L ||
      is.na(num_clusters) || num_clusters < 2) {
    stop("`num_clusters` must be a single integer >= 2.", call. = FALSE)
  }

  structure(
    list(
      num_clusters = as.integer(num_clusters),
      engine = "manydist"
    ),
    class = "pam_dist"
  )
}

#' @exportS3Method base::print
print.pam_dist <- function(x, ...) {
  cat("PAM distance clustering specification\n")
  cat("  num_clusters :", x$num_clusters, "\n")
  cat("  engine       :", x$engine, "\n")
  invisible(x)
}


# fit ----------------------------------------------------------------------

#' @exportS3Method generics::fit
fit.pam_dist <- function(object, x, y = NULL, recipe = NULL, data = NULL, ...) {
  if (!is.null(recipe) && !is.null(data)) {
    return(.fit_pam_dist(spec = object, recipe = recipe, data = data))
  }

  if (inherits(x, "recipe")) {
    if (is.null(data)) {
      stop("When `x` is a recipe, `data` must be supplied.", call. = FALSE)
    }
    return(.fit_pam_dist(spec = object, recipe = x, data = data))
  }

  stop(
    paste0(
      "`fit.pam_dist()` expects either:\n",
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
  cat("  n_obs        :", nrow(x$train_data), "\n")
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

  D_test_medoids <- D_test_train[, object$medoid_idx, drop = FALSE]
  colnames(D_test_medoids) <- paste0("medoid_", seq_len(ncol(D_test_medoids)))

  if (identical(type, "dist")) {
    return(tibble::as_tibble(D_test_medoids))
  }

  assigned <- apply(D_test_medoids, 1, which.min)

  tibble::tibble(.pred_cluster = factor(assigned, levels = seq_len(object$num_clusters)))
}
