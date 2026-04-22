#' Create a LOVO method specification for `compare_lovo_mdist()`
#'
#' Helper to build one method specification for [compare_lovo_mdist()].
#' This allows tidy-style specification of `response`, e.g.
#' `response = Name`, while storing the method definition as a regular list.
#'
#' @param response Optional response column, supplied either unquoted
#'   (e.g. `Name`) or quoted (e.g. `"Name"`).
#' @param ... Additional arguments passed on to [lovo_mdist()] through
#'   [compare_lovo_mdist()].
#'
#' @return A named list of arguments suitable for one element of the
#'   `methods` argument in [compare_lovo_mdist()].
#'
#' @examples
#' \dontrun{
#' methods <- list(
#'   tvd_sup = lovo_method_spec(
#'     response = species,
#'     distance_cat = "tvd",
#'     response_used = TRUE
#'   ),
#'   tvd_unsup = lovo_method_spec(
#'     response = species,
#'     distance_cat = "tvd",
#'     response_used = FALSE
#'   )
#' )
#' }
#'
#' @export
lovo_method_spec <- function(response = NULL, ...) {
  response_name <- NULL
  response_quo <- rlang::enquo(response)

  if (!rlang::quo_is_null(response_quo)) {
    response_expr <- rlang::quo_get_expr(response_quo)

    if (rlang::is_string(response_expr)) {
      response_name <- response_expr
    } else {
      response_name <- rlang::as_name(response_expr)
    }
  }

  c(list(response = response_name), rlang::list2(...))
}
