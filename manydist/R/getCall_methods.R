#' @importFrom stats getCall
NULL

#' @export
#' @method getCall knn_dist
getCall.knn_dist <- function(x, ...) {
  if (!is.null(x$call)) return(x$call)

  rlang::call2(
    "fit_knn_dist",
    x = quote(x),
    y = quote(y),
    k = x$k
  )
}

#' @export
#' @method getCall model_spec
getCall.model_spec <- function(x, ...) {
  # Only handle your custom model; otherwise fall back
  if (!identical(x$method$pkg, "manydist") ||
      !identical(x$method$fun, "nearest_neighbor_dist")) {
    return(NextMethod())
  }

  arg_exprs <- lapply(x$args, rlang::quo_get_expr)

  rlang::call2(
    "nearest_neighbor_dist",
    !!!arg_exprs,
    mode = x$mode
  )
}
