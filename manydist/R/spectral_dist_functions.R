# R/spectral_dist_functions.R

# internal helpers ---------------------------------------------------------

.coerce_square_dist_matrix <- function(D) {
  if (inherits(D, "MDist")) {
    D <- D$distance
  }

  D <- as.matrix(D)

  if (!is.numeric(D)) {
    stop("`D` must be numeric or coercible to a numeric matrix.", call. = FALSE)
  }
  if (nrow(D) != ncol(D)) {
    stop("`D` must be square.", call. = FALSE)
  }
  if (any(!is.finite(D))) {
    stop("`D` must contain only finite values.", call. = FALSE)
  }
  if (max(abs(D - t(D))) > 1e-8) {
    stop("`D` must be symmetric.", call. = FALSE)
  }

  D
}

.affinity_from_dist <- function(D, sigma = NULL) {
  D <- .coerce_square_dist_matrix(D)

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

.fit_spectral_from_recipe <- function(spec, recipe, data) {
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
  D_train <- .coerce_square_dist_matrix(D_train)

  k_val <- rlang::eval_tidy(spec$num_clusters)
  sigma_val <- if (rlang::quo_is_null(spec$sigma)) NULL else rlang::eval_tidy(spec$sigma)
  nstart_val <- rlang::eval_tidy(spec$nstart)

  if (!is.numeric(k_val) || length(k_val) != 1L || is.na(k_val) || k_val < 2) {
    stop("`num_clusters` must evaluate to a single integer >= 2.", call. = FALSE)
  }
  k_val <- as.integer(k_val)

  if (k_val > nrow(D_train)) {
    stop("`num_clusters` cannot exceed the number of training observations.", call. = FALSE)
  }

  if (!is.null(sigma_val) &&
      (!is.numeric(sigma_val) || length(sigma_val) != 1L || is.na(sigma_val) || sigma_val <= 0)) {
    stop("`sigma` must evaluate to NULL or a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(nstart_val) || length(nstart_val) != 1L || is.na(nstart_val) || nstart_val < 1) {
    stop("`nstart` must evaluate to a single integer >= 1.", call. = FALSE)
  }
  nstart_val <- as.integer(nstart_val)

  aff <- .affinity_from_dist(D_train, sigma = sigma_val)
  A <- aff$A
  sigma_val <- aff$sigma

  d <- rowSums(A)
  Dhalf_inv <- diag(ifelse(d > 0, 1 / sqrt(d), 0), nrow(A), ncol(A))

  S <- Dhalf_inv %*% A %*% Dhalf_inv
  ev <- eigen(S, symmetric = TRUE)

  U <- ev$vectors[, seq_len(k_val), drop = FALSE]
  lambda <- ev$values[seq_len(k_val)]
  U_norm <- .row_normalize(U)

  km <- stats::kmeans(U_norm, centers = k_val, nstart = nstart_val)

  structure(
    list(
      spec           = spec,
      fit            = km,
      recipe         = recipe,
      prepped_recipe = prepped_recipe,
      mdist_step     = mdist_step,
      train_matrix   = D_train,
      sigma          = sigma_val,
      train_affinity = A,
      train_degree   = d,
      eigenvalues    = lambda,
      eigenvectors   = U,
      embedding      = U_norm,
      centers        = km$centers,
      clusters       = km$cluster,
      num_clusters   = k_val,
      nstart         = nstart_val
    ),
    class = "spectral_dist_fit"
  )
}

.fit_spectral_from_dist <- function(spec, x) {
  D_train <- .coerce_square_dist_matrix(x)

  k_val <- rlang::eval_tidy(spec$num_clusters)
  sigma_val <- if (rlang::quo_is_null(spec$sigma)) NULL else rlang::eval_tidy(spec$sigma)
  nstart_val <- rlang::eval_tidy(spec$nstart)

  if (!is.numeric(k_val) || length(k_val) != 1L || is.na(k_val) || k_val < 2) {
    stop("`num_clusters` must evaluate to a single integer >= 2.", call. = FALSE)
  }
  k_val <- as.integer(k_val)

  if (k_val > nrow(D_train)) {
    stop("`num_clusters` cannot exceed the number of observations.", call. = FALSE)
  }

  if (!is.null(sigma_val) &&
      (!is.numeric(sigma_val) || length(sigma_val) != 1L || is.na(sigma_val) || sigma_val <= 0)) {
    stop("`sigma` must evaluate to NULL or a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(nstart_val) || length(nstart_val) != 1L || is.na(nstart_val) || nstart_val < 1) {
    stop("`nstart` must evaluate to a single integer >= 1.", call. = FALSE)
  }
  nstart_val <- as.integer(nstart_val)

  aff <- .affinity_from_dist(D_train, sigma = sigma_val)
  A <- aff$A
  sigma_val <- aff$sigma

  d <- rowSums(A)
  Dhalf_inv <- diag(ifelse(d > 0, 1 / sqrt(d), 0), nrow(A), ncol(A))

  S <- Dhalf_inv %*% A %*% Dhalf_inv
  ev <- eigen(S, symmetric = TRUE)

  U <- ev$vectors[, seq_len(k_val), drop = FALSE]
  lambda <- ev$values[seq_len(k_val)]
  U_norm <- .row_normalize(U)

  km <- stats::kmeans(U_norm, centers = k_val, nstart = nstart_val)

  structure(
    list(
      spec           = spec,
      fit            = km,
      recipe         = NULL,
      prepped_recipe = NULL,
      mdist_step     = NULL,
      train_matrix   = D_train,
      sigma          = sigma_val,
      train_affinity = A,
      train_degree   = d,
      eigenvalues    = lambda,
      eigenvectors   = U,
      embedding      = U_norm,
      centers        = km$centers,
      clusters       = km$cluster,
      num_clusters   = k_val,
      nstart         = nstart_val
    ),
    class = "spectral_dist_fit"
  )
}


# specification ------------------------------------------------------------

#' Spectral clustering specification based on manydist dissimilarities
#'
#' @param num_clusters Number of clusters.
#' @param sigma Optional bandwidth for the Gaussian affinity. If `NULL`,
#'   the median pairwise distance is used.
#' @param nstart Number of random starts for k-means.
#'
#' @return A `spectral_dist_spec` object.
#' @export
#' @examples
#'
#'
#' \dontrun{
#' library(manydist)
#' library(palmerpenguins)
#' library(recipes)
#' library(generics)
#'
#' data <- penguins |>
#'   dplyr::select(-species) |>
#'   tidyr::drop_na()
#'
#' rec <- recipes::recipe(~ ., data = data) |>
#'   step_mdist(all_predictors(), preset = "gower", output = "pairwise")
#'
#' spec <- spectral_dist(num_clusters = 3)
#'
#' fit_obj <- generics::fit(spec, recipe = rec, data = data)
#'
#' print(fit_obj)
#' predict(fit_obj)
#' predict(fit_obj, type = "embed")
#' }

spectral_dist <- function(num_clusters = NULL, sigma = NULL, nstart = 50) {
  structure(
    list(
      num_clusters = rlang::enquo(num_clusters),
      sigma = rlang::enquo(sigma),
      nstart = rlang::enquo(nstart),
      engine = "manydist"
    ),
    class = "spectral_dist_spec"
  )
}

#' @exportS3Method base::print
print.spectral_dist_spec <- function(x, ...) {
  k_val <- tryCatch(rlang::eval_tidy(x$num_clusters), error = function(e) NULL)
  sigma_val <- tryCatch(
    if (rlang::quo_is_null(x$sigma)) NULL else rlang::eval_tidy(x$sigma),
    error = function(e) NULL
  )
  nstart_val <- tryCatch(rlang::eval_tidy(x$nstart), error = function(e) NULL)

  cat("Spectral distance clustering specification\n")
  cat("  num_clusters :", if (is.null(k_val)) "<tune/expr>" else k_val, "\n")
  cat("  sigma        :", if (is.null(sigma_val)) "NULL (median heuristic)" else sigma_val, "\n")
  cat("  nstart       :", if (is.null(nstart_val)) "<tune/expr>" else nstart_val, "\n")
  cat("  engine       :", x$engine, "\n")
  invisible(x)
}

#' @exportS3Method parsnip::set_engine
set_engine.spectral_dist_spec <- function(object, engine, ...) {
  object$engine <- engine
  object
}


# fit ----------------------------------------------------------------------

#' @exportS3Method generics::fit
fit.spectral_dist_spec <- function(object, x, y = NULL, recipe = NULL, data = NULL, ...) {
  engine <- object$engine %||% "manydist"

  if (!identical(engine, "manydist")) {
    stop("Only `engine = 'manydist'` is currently supported.", call. = FALSE)
  }

  if (!is.null(recipe) && !is.null(data)) {
    return(.fit_spectral_from_recipe(
      spec = object,
      recipe = recipe,
      data = data
    ))
  }

  if (inherits(x, "recipe")) {
    if (is.null(data)) {
      stop("When `x` is a recipe, `data` must be supplied.", call. = FALSE)
    }
    return(.fit_spectral_from_recipe(
      spec = object,
      recipe = x,
      data = data
    ))
  }

  if (!missing(x) && !is.null(x)) {
    return(.fit_spectral_from_dist(
      spec = object,
      x = x
    ))
  }

  stop(
    paste0(
      "`generics::fit()` for `spectral_dist` expects either:\n",
      "  - `fit(spec, recipe = <recipe>, data = <data>)`, or\n",
      "  - `fit(spec, x = <recipe>, data = <data>)`, or\n",
      "  - `fit(spec, x = <distance matrix/object>)`."
    ),
    call. = FALSE
  )
}


# fitted object methods ----------------------------------------------------

#' @exportS3Method base::print
print.spectral_dist_fit <- function(x, ...) {
  cat("spectral_dist fit\n")
  cat("  num_clusters :", x$num_clusters, "\n")
  cat("  n_obs        :", nrow(x$train_matrix), "\n")
  cat("  sigma        :", x$sigma, "\n")
  cat("  nstart       :", x$nstart, "\n")
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

  if (is.null(object$mdist_step)) {
    stop(
      "`predict.spectral_dist_fit()` with `new_data` is currently only available ",
      "for recipe-based fits using `step_mdist(output = 'pairwise')`.",
      call. = FALSE
    )
  }

  step_obj <- object$mdist_step
  pred_mat <- tibble::as_tibble(new_data)[, step_obj$columns, drop = FALSE]

  D_test_train <- .apply_mdist(step_obj$preprocessor, new_data = pred_mat)
  D_test_train <- as.matrix(D_test_train)

  U_new <- .spectral_embed_new(object, D_test_train)

  if (identical(type, "embed")) {
    return(tibble::as_tibble(U_new, .name_repair = ~ paste0("dim_", seq_along(.x))))
  }

  .spectral_assign_clusters(
    U_new = U_new,
    centers = object$centers,
    num_clusters = object$num_clusters
  )
}
