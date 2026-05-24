#' Spectral clustering from a distance matrix
#'
#' @param D A distance matrix
#' @param k Number of clusters
#' @param affinity_method Method to build affinity
#'
#' @return A vector of cluster labels
#' @export
spectral_from_dist <- function(D, k, affinity_method = "selftune") {
  # source("R/dist_to_affinity.R")
  # source("R/spectral_cluster.R")

   A <- dist_to_affinity(D, method = affinity_method)
  spectral_cluster(A, k)
}
