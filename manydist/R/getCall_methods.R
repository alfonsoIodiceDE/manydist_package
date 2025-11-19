#' @export
#' @method getCall knn_dist
getCall.knn_dist <- function(x, ...) {
  # If we stored the call in fit_knn_dist(), return that
  if (!is.null(x$call)) {
    return(x$call)
  }
  # Fallback: construct a simple synthetic call
  rlang::call2("fit_knn_dist",
               x = quote(x),
               y = quote(y),
               k = x$k)
}

#' @export
#' @method getCall nearest_neighbor_dist
getCall.nearest_neighbor_dist <- function(x, ...) {
  # x is a model_spec; reconstruct a call from its args + mode
  # x$args are quosures; use rlang::quo_get_expr() to get the raw expressions

  arg_exprs <- lapply(x$args, rlang::quo_get_expr)

  rlang::call2(
    "nearest_neighbor_dist",
    !!!arg_exprs,
    mode = x$mode
  )
}
