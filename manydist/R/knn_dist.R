#' Direct nearest-neighbour prediction
#'
#' Predict outcomes from precomputed test-to-training distances, or compute
#' those distances from tabular training and test predictors before predicting.
#'
#' @param x With `dist_fun = NULL`, a numeric test-to-training distance matrix,
#'   data frame, or `MDist` object. Otherwise, tabular training predictors.
#' @param y Training outcomes, in the same order as the distance columns or
#'   training rows. Use a factor for classification and numeric outcomes for
#'   regression. Missing outcomes are not supported. May be omitted when
#'   `response` selects an outcome column from tabular training data.
#' @param response Optional outcome column name (a string or unquoted column
#'   name) in tabular `x`. Supply either `y` or `response`, not both. The selected
#'   column is excluded from the training and test predictors before distance
#'   computation. This selects the prediction outcome, not the `response`
#'   argument of [mdist()].
#' @param new_data Test predictors when `dist_fun` is supplied. Must be `NULL`
#'   when `x` already contains test-to-training distances.
#' @param k Integer number of neighbours, between one and the number of
#'   training observations.
#' @param dist_fun Optional function accepting training predictors as its first
#'   argument and test predictors as `new_data`. It must return a numeric
#'   test-to-training matrix, data frame, or `MDist` object. For example, [mdist()].
#' @param dist_args Named list of additional arguments passed to `dist_fun`.
#' @param type Prediction type: `"class"` for class labels, `"prob"` for class
#'   probabilities, or `"numeric"` for regression predictions.
#'
#' @details
#' Each distance row corresponds to a test observation and each column to a
#' training observation. Distances must be finite and nonnegative. Classification
#' uses an unweighted majority vote; regression uses the mean of the neighbours'
#' outcomes. Distance ties are resolved in training order and class-vote ties
#' in factor-level order. Class probabilities are neighbour class proportions.
#'
#' No training-to-training matrix is required. With `dist_fun = mdist`, distance
#' construction uses `x` as training data and `new_data` as test data. The
#' function computes distances on each call; it does not return a fitted model.
#' Other distance functions are responsible for their own preprocessing rules.
#' For fitted tidymodels workflows, use [nearest_neighbor_dist()] with
#' [step_mdist()] instead.
#'
#' @return A factor for `type = "class"`, a data frame with one column per
#'   outcome level for `type = "prob"`, or a numeric vector for
#'   `type = "numeric"`. Each result has one row or element per test observation.
#'
#' @examples
#' train <- data.frame(value = c(0, 2, 8, 10))
#' test <- data.frame(value = c(1, 9))
#' labels <- factor(c("low", "low", "high", "high"))
#' d <- mdist(train, new_data = test, preset = "gower")
#' knn_dist(d, labels, k = 1)
#' knn_dist(train, labels, new_data = test, k = 1,
#'          dist_fun = mdist, dist_args = list(preset = "gower"))
#' train$label <- labels
#' knn_dist(train, response = "label", new_data = test, k = 1,
#'          dist_fun = mdist, dist_args = list(preset = "gower"))
#' knn_dist(d, labels, k = 2, type = "prob")
#' knn_dist(d, c(0, 2, 8, 10), k = 2, type = "numeric")
#' @seealso [mdist()], [nearest_neighbor_dist()], [step_mdist()]
#' @md
#' @export
knn_dist <- function(x, y = NULL, new_data = NULL, k = 5,
                     dist_fun = NULL, dist_args = list(),
                     type = c("class", "prob", "numeric"), response = NULL) {
  type <- match.arg(type)
  response_quo <- rlang::enquo(response)
  if (!rlang::quo_is_null(response_quo)) {
    if (!is.null(y)) {
      stop("Supply either `y` or `response`, not both.", call. = FALSE)
    }
    if (is.null(dist_fun) || !is.data.frame(x)) {
      stop("`response` requires tabular training data in `x` and a distance function.",
           call. = FALSE)
    }
    response_index <- tidyselect::eval_select(response_quo, x)
    if (length(response_index) != 1L) {
      stop("`response` must select exactly one outcome column.", call. = FALSE)
    }
    response_name <- names(x)[unname(response_index)]
    y <- x[[response_name]]
    x <- x[, setdiff(names(x), response_name), drop = FALSE]
    if (is.data.frame(new_data) && response_name %in% names(new_data)) {
      new_data <- new_data[, setdiff(names(new_data), response_name), drop = FALSE]
    }
  }
  if (!is.null(dim(y)) || length(y) == 0L || anyNA(y)) {
    stop("`y` must be a nonempty vector of nonmissing training outcomes.",
         call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 1 || k != floor(k) || k > length(y)) {
    stop("`k` must be an integer between 1 and the number of training outcomes.",
         call. = FALSE)
  }
  if (type == "numeric") {
    if (!is.numeric(y) || any(!is.finite(y))) {
      stop("Numeric predictions require finite numeric outcomes in `y`.",
           call. = FALSE)
    }
  } else if (!is.factor(y)) {
    stop("Classification requires factor outcomes in `y`.", call. = FALSE)
  }
  if (!is.list(dist_args) ||
      (length(dist_args) > 0L &&
       (is.null(names(dist_args)) || anyNA(names(dist_args)) ||
        any(names(dist_args) == "") || anyDuplicated(names(dist_args))))) {
    stop("`dist_args` must be a named list with unique, nonempty names.",
         call. = FALSE)
  }
  if (any(c("x", "new_data") %in% names(dist_args))) {
    stop("Supply training and test data through `x` and `new_data`, not `dist_args`.",
         call. = FALSE)
  }

  expected_rows <- NULL
  if (is.null(dist_fun)) {
    if (!is.null(new_data) || length(dist_args) > 0L) {
      stop("With `dist_fun = NULL`, supply distances in `x` and leave `new_data` and `dist_args` empty.",
           call. = FALSE)
    }
    distances <- x
  } else {
    if (!is.function(dist_fun)) {
      stop("`dist_fun` must be a function or NULL.", call. = FALSE)
    }
    if (!(is.data.frame(x) || is.matrix(x)) || nrow(x) != length(y)) {
      stop("Training predictors in `x` must have one row per outcome in `y`.",
           call. = FALSE)
    }
    if (!(is.data.frame(new_data) || is.matrix(new_data))) {
      stop("Supply tabular test predictors in `new_data` when using `dist_fun`.",
           call. = FALSE)
    }
    expected_rows <- nrow(new_data)
    distances <- do.call(dist_fun, c(list(x, new_data = new_data), dist_args))
  }

  if (inherits(distances, "MDist")) distances <- distances$distance
  if (!(is.matrix(distances) || is.data.frame(distances))) {
    stop("Distances must be a numeric test-to-training matrix, data frame, or MDist object.",
         call. = FALSE)
  }
  D <- as.matrix(distances)
  if (!is.numeric(D) || any(!is.finite(D)) || any(D < 0)) {
    stop("Distances must be numeric, finite, and nonnegative.", call. = FALSE)
  }
  if (ncol(D) != length(y)) {
    stop("The distance matrix must have one column per training outcome in `y`.",
         call. = FALSE)
  }
  if (!is.null(expected_rows) && nrow(D) != expected_rows) {
    stop("The distance function must return one row per test observation.",
         call. = FALSE)
  }
  k <- as.integer(k)
  switch(type,
         class = .knn_class_from_dist(D, y, k),
         prob = .knn_prob_from_dist(D, y, k),
         numeric = .knn_reg_from_dist(D, y, k))
}
