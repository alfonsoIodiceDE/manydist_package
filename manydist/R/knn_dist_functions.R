# Internal helpers (NOT exported) ------------------------------

# majority vote (classification)
.knn_class_from_dist <- function(D, y_train, k) {
  nn_idx <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  lv <- levels(y_train)
  out <- apply(nn_idx, 2L, function(idx) {
    tab <- table(y_train[idx])
    names(tab)[which.max(tab)]
  })
  factor(out, levels = lv)
}

# posterior probabilities
.knn_prob_from_dist <- function(D, y_train, k) {
  classes <- levels(y_train)
  nn_idx <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  probs <- t(apply(nn_idx, 2L, function(idx) {
    tab <- table(factor(y_train[idx], levels = classes))
    as.numeric(tab / sum(tab))
  }))
  colnames(probs) <- classes
  as.data.frame(probs)
}

# regression mean
.knn_reg_from_dist <- function(D, y_train, k) {
  nn_idx <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  apply(nn_idx, 2L, function(idx) mean(y_train[idx]))
}

#' KNN fit with precomputed distances
#'
#' This is the engine function used by the parsnip
#' `nearest_neighbor_dist(precomputed)` model.
#'
#' @param x A data frame or matrix of predictors (already a distance
#'   representation if `dist_fun` is NULL).
#' @param y A vector of outcomes.
#' @param k Number of neighbors.
#' @param dist_fun A distance function taking arguments `(x, new_data = ...)`,
#'   or NULL if `x`/`new_data` are already distance matrices.
#' @param dist_args A named list of additional arguments passed to `dist_fun`.
#'
#' @export
fit_knn_dist <- function(x, y, k = 5, dist_fun = NULL, dist_args = list()) {

  if (!is.null(dist_fun) && !is.function(dist_fun)) {
    stop("`dist_fun` must be a function or NULL.")
  }

  if (k < 1) {
    stop("`k` must be >= 1.")
  }

  # store the call for downstream tools that expect a `$call` slot
  cl <- match.call()

  structure(
    list(
      x         = x,
      y         = y,
      k         = k,
      dist_fun  = dist_fun,
      dist_args = dist_args,
      call      = match.call()   # <-- THIS is the crucial addition
    ),
    class = "knn_dist"
  )
}
#' Predict from a KNN-distance model (classification)
#'
#' Used by parsnip for `type = "class"` and `type = "prob"`.
#'
#' @param object A `knn_dist` object created by `fit_knn_dist()`.
#' @param new_data A data frame or matrix of precomputed distances
#'   *or* raw predictors, depending on `dist_fun`.
#' @param type `"class"` or `"prob"`.
#'
#' @export
predict_knn_dist_class <- function(object, new_data, type = c("class", "prob")) {
  type <- match.arg(type)

  if (is.null(object$dist_fun)) {

    # Expected: new_data already IS a distance matrix
    D <- as.matrix(new_data)

  } else {

    D_raw <- do.call(
      object$dist_fun,
      c(list(object$x, new_data = new_data), object$dist_args)
    )

    D <- if (inherits(D_raw, "MDist")) {
      as.matrix(D_raw$distance)
    } else {
      D_raw
    }
  }

  if (type == "class") .knn_class_from_dist(D, object$y, object$k)
  else                 .knn_prob_from_dist(D, object$y, object$k)
}

#' Predict class probabilities from a KNN-distance model
#'
#' Thin wrapper used by parsnip for `type = "prob"`.
#'
#' @param object A `knn_dist` object created by `fit_knn_dist()`.
#' @param new_data A data frame or matrix of precomputed distances
#'   *or* raw predictors, depending on `dist_fun`.
#' @export
predict_knn_dist_prob <- function(object, new_data, ...) {
  predict_knn_dist_class(object, new_data, type = "prob")
}

#' Predict from a KNN-distance model (regression)
#'
#' Used by parsnip for regression mode (`type = "numeric"`).
#' Supports two modes:
#' - `dist_fun` is NULL  →  `new_data` is assumed to already be a
#'   distance matrix (e.g., produced by `step_mdist()`).
#' - `dist_fun` is a function  →  distances are computed as
#'   `dist_fun(x = object$x, new_data = new_data, ...)`.
#'
#' @param object A `knn_dist` object created by `fit_knn_dist()`.
#' @param new_data A data frame or matrix. Either:
#'   - precomputed distances to training points (if `dist_fun` is NULL), or
#'   - raw predictors to be passed to `dist_fun`.
#'
#' @export
predict_knn_dist_reg <- function(object, new_data) {
  # Case 1: precomputed distances → new_data *is* the distance matrix
  if (is.null(object$dist_fun)) {
    D <- as.matrix(new_data)

  } else {
    # Case 2: need to compute distances via dist_fun
    D_raw <- do.call(
      object$dist_fun,
      c(list(object$x, new_data = new_data), object$dist_args)
    )

    if (inherits(D_raw, "MDist")) {
      D <- as.matrix(D_raw$distance)
    } else if (is.matrix(D_raw)) {
      D <- D_raw
    } else {
      stop("dist_fun must return either an 'MDist' object or a numeric matrix.")
    }
  }

  .knn_reg_from_dist(D, object$y, object$k)
}
