.affinity_from_dist <- function(D, sigma = NULL) {
  D <- as.matrix(D)

  if (!is.numeric(D)) {
    stop("D must be numeric.", call. = FALSE)
  }
  if (nrow(D) != ncol(D)) {
    stop("D must be square.", call. = FALSE)
  }
  if (any(!is.finite(D))) {
    stop("D must contain only finite values.", call. = FALSE)
  }
  if (max(abs(D - t(D))) > 1e-8) {
    stop("D must be symmetric.", call. = FALSE)
  }

  D <- (D + t(D)) / 2
  diag(D) <- 0

  if (is.null(sigma)) {
    sigma <- stats::median(D[upper.tri(D)])
  }

  if (!is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive finite value.", call. = FALSE)
  }

  A <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0
  A <- (A + t(A)) / 2

  list(A = A, sigma = sigma)
}

.affinity_new_to_train <- function(D_new, sigma) {
  D_new <- as.matrix(D_new)

  if (!is.numeric(D_new)) {
    stop("`D_new` must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(D_new))) {
    stop("`D_new` must contain only finite values.", call. = FALSE)
  }
  if (!is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive finite value.", call. = FALSE)
  }

  exp(-(D_new^2) / (2 * sigma^2))
}

.row_normalize <- function(U) {
  norms <- sqrt(rowSums(U^2))
  norms[norms == 0] <- 1
  U / norms
}

.spectral_embed_new <- function(object, D_new) {
  D_new <- as.matrix(D_new)

  if (ncol(D_new) != nrow(object$train_matrix)) {
    stop("`D_new` must have one column for each training observation.", call. = FALSE)
  }

  A_new <- .affinity_new_to_train(D_new, sigma = object$sigma)
  d_new <- rowSums(A_new)

  S_new <- sweep(A_new, 1, sqrt(pmax(d_new, 1e-12)), "/")
  S_new <- sweep(S_new, 2, sqrt(pmax(object$train_degree, 1e-12)), "/")

  U_new <- S_new %*% object$eigenvectors
  U_new <- sweep(U_new, 2, pmax(object$eigenvalues, 1e-12), "/")

  .row_normalize(U_new)
}

.spectral_assign_clusters <- function(U_new, centers, num_clusters) {
  assigned <- apply(U_new, 1, function(x) {
    d <- colSums((t(centers) - x)^2)
    which.min(d)
  })

  tibble::tibble(.pred_cluster = factor(assigned, levels = seq_len(num_clusters)))
}

.fit_spectral_dist <- function(spec, recipe, data) {
  if (missing(recipe) || missing(data)) {
    stop("`recipe` and `data` must be provided.", call. = FALSE)
  }

  prepped_recipe <- recipes::prep(recipe, training = data)

  mdist_step <- .get_step_mdist(prepped_recipe)

  if (!identical(mdist_step$output, "pairwise")) {
    stop(
      "`spectral_dist` requires `step_mdist(..., output = 'pairwise')` in the recipe.",
      call. = FALSE
    )
  }

  baked_train <- recipes::bake(prepped_recipe, new_data = data)
  D_train <- .extract_dist_block(baked_train)

  if (nrow(D_train) != ncol(D_train)) {
    stop("The baked training dissimilarity matrix is not square.", call. = FALSE)
  }

  D_train <- unname(D_train)
  D_train <- (D_train + t(D_train)) / 2
  diag(D_train) <- 0

  if (spec$num_clusters < 2 || spec$num_clusters > nrow(D_train)) {
    stop(
      "`num_clusters` must be between 2 and the number of training observations.",
      call. = FALSE
    )
  }

  aff <- .affinity_from_dist(D_train, sigma = spec$sigma)
  A <- aff$A
  sigma <- aff$sigma

  n <- nrow(A)
  d <- rowSums(A)
  Dhalf_inv <- diag(ifelse(d > 0, 1 / sqrt(d), 0), n, n)

  S <- Dhalf_inv %*% A %*% Dhalf_inv
  ev <- eigen(S, symmetric = TRUE)

  U <- ev$vectors[, seq_len(spec$num_clusters), drop = FALSE]
  lambda <- ev$values[seq_len(spec$num_clusters)]
  U_norm <- .row_normalize(U)

  km <- stats::kmeans(U_norm, centers = spec$num_clusters, nstart = spec$nstart)

  structure(
    list(
      spec           = spec,
      recipe         = recipe,
      prepped_recipe = prepped_recipe,
      mdist_step     = mdist_step,
      train_data     = tibble::as_tibble(data),
      train_matrix   = D_train,
      sigma          = sigma,
      train_affinity = A,
      train_degree   = d,
      eigenvalues    = lambda,
      eigenvectors   = U,
      embedding      = U_norm,
      kmeans         = km,
      clusters       = km$cluster,
      num_clusters   = spec$num_clusters
    ),
    class = "spectral_dist_fit"
  )
}


# specification ------------------------------------------------------------

#' Spectral clustering specification for manydist dissimilarities
#'
#' @param num_clusters Number of clusters.
#' @param sigma Optional bandwidth for the Gaussian affinity. If `NULL`,
#'   the median pairwise distance is used.
#' @param nstart Number of random starts for k-means.
#'
#' @return A `spectral_dist` model specification.
#' @export
spectral_dist <- function(num_clusters, sigma = NULL, nstart = 50) {
  if (missing(num_clusters)) {
    stop("`num_clusters` must be provided.", call. = FALSE)
  }

  if (!is.numeric(num_clusters) || length(num_clusters) != 1L ||
      is.na(num_clusters) || num_clusters < 2) {
    stop("`num_clusters` must be a single integer >= 2.", call. = FALSE)
  }

  if (!is.null(sigma) &&
      (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || sigma <= 0)) {
    stop("`sigma` must be `NULL` or a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(nstart) || length(nstart) != 1L || is.na(nstart) || nstart < 1) {
    stop("`nstart` must be a single integer >= 1.", call. = FALSE)
  }

  structure(
    list(
      num_clusters = as.integer(num_clusters),
      sigma = sigma,
      nstart = as.integer(nstart),
      engine = "manydist"
    ),
    class = "spectral_dist"
  )
}

#' @exportS3Method base::print
print.spectral_dist <- function(x, ...) {
  cat("Spectral distance clustering specification\n")
  cat("  num_clusters :", x$num_clusters, "\n")
  cat("  sigma        :", if (is.null(x$sigma)) "NULL (median heuristic)" else x$sigma, "\n")
  cat("  nstart       :", x$nstart, "\n")
  cat("  engine       :", x$engine, "\n")
  invisible(x)
}


# fit ----------------------------------------------------------------------

#' @exportS3Method generics::fit
fit.spectral_dist <- function(object, x, y = NULL, recipe = NULL, data = NULL, ...) {
  if (!is.null(recipe) && !is.null(data)) {
    return(.fit_spectral_dist(spec = object, recipe = recipe, data = data))
  }

  if (inherits(x, "recipe")) {
    if (is.null(data)) {
      stop("When `x` is a recipe, `data` must be supplied.", call. = FALSE)
    }
    return(.fit_spectral_dist(spec = object, recipe = x, data = data))
  }

  stop(
    paste0(
      "`fit.spectral_dist()` expects either:\n",
      "  - `fit(spec, recipe = <recipe>, data = <data>)`, or\n",
      "  - `fit(spec, x = <recipe>, data = <data>)`."
    ),
    call. = FALSE
  )
}


# fitted object methods ----------------------------------------------------

#' @exportS3Method base::print
print.spectral_dist_fit <- function(x, ...) {
  cat("spectral_dist fit\n")
  cat("  num_clusters :", x$num_clusters, "\n")
  cat("  n_obs        :", nrow(x$train_data), "\n")
  cat("  sigma        :", x$sigma, "\n")
  invisible(x)
}

#' @exportS3Method stats::predict
predict.spectral_dist_fit <- function(object, new_data = NULL, type = c("cluster", "embed"), ...) {
  type <- match.arg(type)

  if (is.null(new_data)) {
    if (identical(type, "cluster")) {
      return(tibble::tibble(.pred_cluster = factor(object$clusters)))
    }
    return(tibble::as_tibble(object$embedding, .name_repair = ~ paste0("dim_", seq_along(.x))))
  }

  step_obj <- object$mdist_step
  new_tbl <- tibble::as_tibble(new_data)

  if (!all(step_obj$columns %in% names(new_tbl))) {
    stop("`new_data` must contain the same predictor columns used by `step_mdist()`.", call. = FALSE)
  }

  pred_mat <- new_tbl[, step_obj$columns, drop = FALSE]

  D_test_train <- .apply_mdist(step_obj$preprocessor, new_data = pred_mat)
  D_test_train <- as.matrix(D_test_train)

  U_new <- .spectral_embed_new(object, D_test_train)

  if (identical(type, "embed")) {
    return(tibble::as_tibble(U_new, .name_repair = ~ paste0("dim_", seq_along(.x))))
  }

  .spectral_assign_clusters(
    U_new = U_new,
    centers = object$kmeans$centers,
    num_clusters = object$num_clusters
  )
}
