gudmm_jensen_shannon <- function(p, q, base = NULL) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  # Normalize if they don't sum to 1
  if (sum(p) > 0) p <- p / sum(p)
  if (sum(q) > 0) q <- q / sum(q)
  
  m <- (p + q) / 2.0
  
  # Compute relative entropy (KL divergence)
  left <- ifelse(p == 0 | m == 0, 0, p * log(p / m))
  right <- ifelse(q == 0 | m == 0, 0, q * log(q / m))
  
  js <- sum(left) + sum(right)
  
  if (!is.null(base)) {
    js <- js / log(base)
  }
  
  return(sqrt(js / 2.0))
}
