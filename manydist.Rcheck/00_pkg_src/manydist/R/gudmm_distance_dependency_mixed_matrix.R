# Main distance calculation function (exact translation)
gudmm_distance_dependency_mixed_matrix <- function(XX, no_f_cont, no_f_ord, DM) {
    
  
  X <- XX
  N <- nrow(X)
  no_f <- ncol(X)
  
  # Remove features with no variation - exact as Python
  rep_fea <- c()
  for (i in 1:no_f) {
    unique_vals <- unique(X[, i][!is.na(X[, i])])
    if (length(unique_vals) == 1) {
      rep_fea <- c(rep_fea, i)
    }
  }
  
  if (length(rep_fea) > 0) {
    X <- X[, -rep_fea, drop = FALSE]
    no_f <- no_f - length(rep_fea)
    no_f_ord <- no_f_ord - length(rep_fea)
  }
  
  no_f_disc <- no_f - no_f_cont
  
  # Scale continuous features - exact as Python (MinMaxScaler equivalent)
  if (no_f_cont > 0) {
    for (i in 1:no_f_cont) {
      col_data <- X[, i]
      non_na_data <- col_data[!is.na(col_data)]
      if (length(non_na_data) > 1) {
        min_val <- min(non_na_data)
        max_val <- max(non_na_data)
        if (max_val != min_val) {
          X[!is.na(col_data), i] <- (non_na_data - min_val) / (max_val - min_val)
        }
      }
    }
  }
  
 
  # Calculate mutual information matrix - exact
  R <- gudmm_calculate_R_MI(X, no_f_cont)

  # Calculate categorical distances - exact
  if (DM == 'DM0') {
    D_rth <- NULL
  } else {
    D_rth <- gudmm_distance_categorical_exact(X, no_f_cont, R, DM)
  }
  
  # Calculate pairwise distances for discrete part - exact translation
  dist_discrete <- matrix(0, N, N)
  
  for (i in 1:N) {
    for (j in i:N) {
      if (DM == 'DM0') {
        # Hamming distance - exact as Python
        non_na_mask <- !is.na(X[i, ]) & !is.na(X[j, ])
        if (sum(non_na_mask) > 0) {
          hamming_dist <- sum(X[i, non_na_mask] != X[j, non_na_mask]) / sum(non_na_mask)
          dist_discrete[i, j] <- (no_f_disc / (no_f_disc + 1)) * hamming_dist
        }
      } else {
        if (no_f_disc > 0 && !is.null(D_rth) && length(D_rth) > 0) {
          dist_ij_r <- 0  # Initialize as scalar
          
          for (r in (no_f_cont + 1):no_f) {
            r_idx <- r - no_f_cont
            
            if (r_idx <= length(D_rth) && !is.na(X[i, r]) && !is.na(X[j, r])) {
              x_i <- as.integer(X[i, r])
              x_j <- as.integer(X[j, r])
              
              # Ensure valid indices (1-based in R)
              if (x_i > 0 && x_j > 0 && 
                  x_i <= nrow(D_rth[[r_idx]]) && 
                  x_j <= ncol(D_rth[[r_idx]])) {
                dist_ij_r <- dist_ij_r + D_rth[[r_idx]][x_i, x_j]
              }
            }
          }
          
          dist_discrete[i, j] <- dist_ij_r
        }
      }
      dist_discrete[j, i] <- dist_discrete[i, j]
    }
  }
  
  # Calculate continuous distances - exact translation
  if (no_f_cont > 0) {
    cont_data <- X[, 1:no_f_cont, drop = FALSE]
    
    if (DM == 'DM4') {
      # Original Mahalanobis - exact as Python
      if (no_f_cont == 1) {
        VI <- 1 / var(cont_data, na.rm = TRUE)
      } else {
        cov_mat <- cov(cont_data, use = "complete.obs")
        VI <- solve(cov_mat)
      }
      
      # Calculate Mahalanobis distances
      MD <- matrix(0, N, N)
      for (i in 1:N) {
        for (j in 1:N) {
          diff_vec <- cont_data[i, ] - cont_data[j, ]
          if (!any(is.na(diff_vec))) {
            if (no_f_cont == 1) {
              MD[i, j] <- sqrt(abs(diff_vec^2 * VI))
            } else {
              MD[i, j] <- sqrt(abs(t(diff_vec) %*% VI %*% diff_vec))
            }
          } else {
            MD[i, j] <- sqrt(sum(diff_vec^2, na.rm = TRUE))
          }
        }
      }
      
    } else if (DM == 'DM5') {
      # Modified Mahalanobis with R weighting - exact as Python
      if (no_f_cont == 1) {
        var_val <- var(cont_data, na.rm = TRUE)
        if (var_val > 0) {
          VI <- (1 / var_val) * R[1:no_f_cont, 1:no_f_cont]
        } else {
          VI <- R[1:no_f_cont, 1:no_f_cont]
        }
      } else {
        cov_mat <- cov(cont_data, use = "complete.obs")
        if (det(cov_mat) > 1e-10) {
          VI <- solve(cov_mat) * R[1:no_f_cont, 1:no_f_cont]
        } else {
          VI <- diag(no_f_cont) * diag(R[1:no_f_cont, 1:no_f_cont])
        }
      }
      
      # Calculate weighted Mahalanobis distances
      MD <- matrix(0, N, N)
      for (i in 1:N) {
        for (j in 1:N) {
          diff_vec <- cont_data[i, ] - cont_data[j, ]
          if (!any(is.na(diff_vec))) {
            if (no_f_cont == 1) {
              MD[i, j] <- sqrt(abs(diff_vec^2 * VI))
            } else {
              MD[i, j] <- sqrt(abs(t(diff_vec) %*% VI %*% diff_vec))
            }
          } else {
            MD[i, j] <- sqrt(sum(diff_vec^2, na.rm = TRUE))
          }
        }
      }
      
    } else {
      # Euclidean distance - exact as Python
      MD <- as.matrix(dist(cont_data, method = "euclidean"))
    }
  } else {
    MD <- matrix(0, N, N)
  }
  
  dist_continuous <- MD
  
  # Combine distances - exact formula from Python
  final_dist <- (1 / (no_f_disc + 1)) * dist_continuous + 
    (no_f_disc / (no_f_disc + 1)) * dist_discrete
  
  return(final_dist)
}
