dist_to_affinity <- function(D, method = c("gaussian", "selftune"), sigma = NULL, k = 7) {
  method <- match.arg(method)
  n <- nrow(D)
  
  if (method == "gaussian") {
    if (is.null(sigma)) sigma <- median(as.dist(D))
    A <- exp(-(D^2) / (2 * sigma^2))
    diag(A) <- 0
    return(A)
  } else {
    # self-tuning affinity (Zelnik–Manor & Perona)
    kth <- apply(D, 1, function(row) sort(row[row > 0])[k])
    A <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) A[i,j] <- exp(-(D[i,j]^2) / (kth[i] * kth[j]))
      }
    }
    return((A + t(A)) / 2)
  }
}
