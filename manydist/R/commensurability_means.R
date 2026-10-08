# Mean over the full training square, including the diagonal, matching the
# existing train-train convention. No observation-by-observation matrix needed.
.mean_absolute_training_distance <- function(x) {
  x <- sort(as.numeric(x), na.last = TRUE)
  n <- length(x)
  if (n == 0L || any(!is.finite(x))) return(NA_real_)
  x <- x - x[1L]
  2 * sum((2 * seq_len(n) - n - 1) * x) / n^2
}

.mean_categorical_training_distance <- function(Z, delta) {
  frequencies <- colMeans(as.matrix(Z))
  as.numeric(crossprod(frequencies, as.matrix(delta) %*% frequencies))
}

.mean_indicator_training_distance <- function(Z, eta) {
  prototypes <- diag(eta, nrow = length(eta))
  delta <- as.matrix(stats::dist(prototypes, method = "euclidean"))
  .mean_categorical_training_distance(Z, delta)
}

.commensurability_denominator <- function(mean_distance) {
  if (!is.finite(mean_distance) || mean_distance < 0) {
    stop("Cannot estimate a finite nonnegative training mean distance.",
         call. = FALSE)
  }
  # An identically zero training contribution has no estimable positive scale.
  if (mean_distance == 0) 1 else mean_distance
}
