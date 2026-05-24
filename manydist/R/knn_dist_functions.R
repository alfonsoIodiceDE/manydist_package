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

#' @rdname nearest_neighbor_dist
#' @keywords internal
#'
#' @param x A data frame or matrix of predictors. If `dist_fun = NULL`,
#'   `x` is treated as an already computed distance representation.
#' @param y A vector of outcomes.
#' @param k Number of neighbours.
#' @param dist_fun Optional distance function taking arguments `x` and
#'   `new_data`. If `NULL`, inputs are assumed to already be distance matrices.
#' @param dist_args A named list of additional arguments passed to `dist_fun`.
#'
#' @export
fit_knn_dist <- function(x, y, k = 5, dist_fun = NULL, dist_args = list()) {

  if (rlang::is_quosure(k) || inherits(k, "formula")) {
    k <- rlang::eval_tidy(k)
  }

  k <- as.integer(k)

  if (!is.null(dist_fun) && !is.function(dist_fun)) {
    stop("`dist_fun` must be a function or NULL.")
  }

  if (length(k) != 1 || is.na(k) || k < 1) {
    stop("`k` must be a single integer >= 1.")
  }

  structure(
    list(
      x         = x,
      y         = y,
      k         = k,
      dist_fun  = dist_fun,
      dist_args = dist_args,
      call      = match.call()
    ),
    class = "knn_dist"
  )
}

#' @rdname nearest_neighbor_dist
#' @keywords internal
#'
#' @param object A fitted `knn_dist` object created by `fit_knn_dist()`.
#' @param new_data A data frame or matrix of precomputed distances, or raw
#'   predictors when `dist_fun` is supplied.
#' @param type Prediction type. For classification, `"class"` returns class
#'   predictions and `"prob"` returns class probabilities.
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

#' @rdname nearest_neighbor_dist
#' @keywords internal
#'
#' @param ... Additional arguments passed to lower-level methods or currently
#'   not used.
#'
#' @export
predict_knn_dist_prob <- function(object, new_data, ...) {
  predict_knn_dist_class(object, new_data, type = "prob")
}

#' @rdname nearest_neighbor_dist
#' @keywords internal
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
