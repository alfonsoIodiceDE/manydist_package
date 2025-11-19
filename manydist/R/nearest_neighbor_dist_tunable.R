#' Tuning parameters for nearest_neighbor_dist
#'
#' This defines which arguments of `nearest_neighbor_dist()` are tunable
#' and which dials parameter to use for each one.
#'
#' @param x A `nearest_neighbor_dist` model specification.
#' @param ... Not used (included for method consistency).
#'
#' @export
#' @method tunable nearest_neighbor_dist
tunable.nearest_neighbor_dist <- function(x, ...) {
  tibble::tibble(
    name = "neighbors",

    # dials param object (this is what hardhat/parameters() really want)
    object = list(dials::neighbors()),

    # optional; still useful for documentation
    call_info = list(
      list(
        pkg = "dials",
        fun = "neighbors"
      )
    ),

    source       = "model_spec",
    component    = "nearest_neighbor_dist",
    component_id = "main"
  )
}
