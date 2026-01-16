.onLoad <- function(libname, pkgname) {
  if (requireNamespace("parsnip", quietly = TRUE)) {
    try(register_nearest_neighbor_dist(), silent = TRUE)
  }
}
