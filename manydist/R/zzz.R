.onLoad <- function(libname, pkgname) {
  if (requireNamespace("parsnip", quietly = TRUE)) {
    try(manydist:::register_nearest_neighbor_dist(), silent = TRUE)
  }
}
