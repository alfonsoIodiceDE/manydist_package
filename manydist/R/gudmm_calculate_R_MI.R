# Calculate mutual information matrix R (exact translation)
gudmm_calculate_R_MI <- function(X, no_f_cont) {
  source("R/gudmm_mutual_info_regression.R")
  source("R/gudmm_mutual_info_classif.R")
  source("R/gudmm_mutual_info_score.R")
  
  no_f <- ncol(X)
  R <- matrix(0, no_f, no_f)
  
  for (r in 1:no_f) {
    x_r <- X[, r, drop = FALSE]
    for (s in r:no_f) {
      x_s <- X[, s, drop = FALSE]
      
      if (s <= no_f_cont && r <= no_f_cont) {
        # Both continuous - mutual information regression
        R[r, s] <- gudmm_mutual_info_regression(x_r, x_s)
      } else if (r <= no_f_cont && s > no_f_cont) {
        # r continuous, s categorical - mutual information classification
        R[r, s] <- gudmm_mutual_info_classif(x_r, x_s)
      } else if (r > no_f_cont && s > no_f_cont) {
        # Both categorical - mutual information score
        R[r, s] <- gudmm_mutual_info_score(x_r, x_s)
      }
      
      R[s, r] <- R[r, s]
    }
  }
  
  # Normalize by diagonal
  diag_r <- diag(R)
  for (r in 1:no_f) {
    if (diag_r[r] > 0) {
      R[r, ] <- R[r, ] / diag_r[r]
    }
  }
  
  return(R)
}

