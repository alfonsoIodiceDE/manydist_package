spectral_cluster <- function(A, k) {
  n <- nrow(A)
  d <- rowSums(A)
  Dhalf_inv <- diag(ifelse(d > 0, 1 / sqrt(d), 0))
  L <- diag(n) - Dhalf_inv %*% A %*% Dhalf_inv
  ev <- eigen(L, symmetric = TRUE)
  U <- ev$vectors[, (n - k + 1):n, drop = FALSE]
  U <- U / sqrt(rowSums(U^2))
  kmeans(U, centers = k, nstart = 50)$cluster
}
