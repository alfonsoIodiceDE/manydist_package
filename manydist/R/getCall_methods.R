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
#' @method getCall nearest_neighbor_dist
getCall.nearest_neighbor_dist <- function(x, ...) {
  arg_exprs <- lapply(x$args, rlang::quo_get_expr)

  rlang::call2(
    "nearest_neighbor_dist",
    !!!arg_exprs,
    mode = x$mode
  )
}
