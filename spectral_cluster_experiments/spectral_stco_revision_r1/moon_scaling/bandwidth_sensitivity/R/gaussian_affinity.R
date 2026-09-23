# Global Gaussian affinity with an explicit multiplier on a reproducible base
# scale. Source the frozen full-study spectral engine before this file.

moon_bandwidth_gaussian_affinity <- function(distance, multiplier = 1) {
  if (!exists("full_validate_distance", mode = "function")) {
    stop("Source the frozen spectral engine first.", call. = FALSE)
  }
  distance <- full_validate_distance(distance)
  if (length(multiplier) != 1L || !is.finite(multiplier) || multiplier <= 0) {
    stop("`multiplier` must be one finite positive number.", call. = FALSE)
  }
  positive <- distance[upper.tri(distance) & distance > 0]
  base_scale <- stats::median(positive)
  scale <- multiplier * base_scale
  if (!length(positive) || !is.finite(scale) || scale <= 0) {
    stop("Gaussian affinity has no positive finite scale.", call. = FALSE)
  }
  affinity <- exp(-(distance^2) / (2 * scale^2))
  diag(affinity) <- 0
  list(
    matrix = affinity,
    base_scale = base_scale,
    scale = scale,
    multiplier = multiplier
  )
}

moon_bandwidth_score_distance <- function(
    distance,
    truth,
    multiplier,
    cluster_seed,
    kmeans_specification
) {
  affinity <- moon_bandwidth_gaussian_affinity(distance, multiplier)
  fit <- full_spectral_cluster(
    affinity$matrix,
    k = nlevels(factor(truth)),
    seed = cluster_seed,
    kmeans_specification = kmeans_specification
  )
  list(
    ARI = mclust::adjustedRandIndex(fit$cluster, truth),
    objective = fit$objective,
    iterations = fit$iterations,
    ifault = fit$ifault,
    affinity_base_scale = affinity$base_scale,
    affinity_scale = affinity$scale,
    affinity_multiplier = affinity$multiplier
  )
}

