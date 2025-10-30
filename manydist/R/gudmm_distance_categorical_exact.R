gudmm_distance_categorical_exact <- function(X, no_f_cont, R, DM) {
  
  
  
  N <- nrow(X)
  no_f <- ncol(X)
  n_sampling <- 1001
  no_f_disc <- no_f - no_f_cont
  
  if (no_f_disc <= 0) {
    return(list())
  }
  
  # Initialization - exact translation
  max_cats <- numeric(no_f_disc)
  phi_rth <- list()
  
  for (r in 1:no_f_disc) {
    max_cats[r] <- max(X[, r + no_f_cont], na.rm = TRUE)
    phi_rth[[r]] <- matrix(0, max_cats[r], max_cats[r])
  }
  
  dist_r_s_th <- list()
  dist_rth <- list()
  D_rth <- list()
  
  for (r in 1:no_f_disc) {
    max_cat_r <- max_cats[r]
    dist_r_s_th[[r]] <- list()
    dist_rth[[r]] <- matrix(0, max_cat_r, max_cat_r)
    D_rth[[r]] <- matrix(0, max_cat_r, max_cat_r)
    
    for (s in 1:no_f_cont) {
      dist_r_s_th[[r]][[s]] <- matrix(0, max_cat_r, max_cat_r)
    }
  }
  
  # Distance Calculation - exact translation
  S_As <- log2(max_cats)
  
  for (r in 1:no_f) {
    if (r > no_f_cont) {
      max_cat_r <- max(X[, r], na.rm = TRUE)
      
      if (DM != 'DM1' && DM != 'DM2') {
        # Continuous-Categorical interaction distances using KDE
        for (s in 1:no_f_cont) {
          if (R[r, s] > 0) {
            KDEs <- matrix(0, n_sampling, max_cat_r)
            x_s <- X[, s]
            
            # Remove NA values for range calculation
            x_s_clean <- x_s[!is.na(x_s)]
            if (length(x_s_clean) == 0) next
            
            a1 <- min(x_s_clean)
            a2 <- max(x_s_clean)
            sm <- sd(x_s_clean)
            ni <- length(x_s_clean)
            
            if (is.na(sm) || sm == 0) sm <- 0.01
            if (a1 == a2) {
              a1 <- a1 - 0.1
              a2 <- a2 + 0.1
            }
            
            h0 <- 1.06 * sm * ni^(-1/5)
            X_plot <- seq(a1, a2, length.out = n_sampling)
            
            # For each category value, estimate KDE
            for (i in 1:max_cat_r) {
              logic_ri <- (X[, r] == i) & !is.na(X[, r]) & !is.na(X[, s])
              p_r_i <- sum(logic_ri) / N
              
              if (sum(logic_ri) > 0) {
                X_s_ri <- X[logic_ri, s] * p_r_i  # Exact as in Python
                X_s_ri_clean <- X_s_ri[!is.na(X_s_ri)]
                
                if (length(X_s_ri_clean) > 1) {
                  sm_sr <- sd(X_s_ri_clean)
                  ni <- length(X_s_ri_clean)
                  
                  if (is.na(sm_sr) || sm_sr == 0) {
                    sm_sr <- 0.01
                  }
                  
                  h0 <- 1.06 * sm_sr * ni^(-1/5)
                  
                  # Kernel Density Estimation - exact translation
                  tryCatch({
                    kde_result <- density(X_s_ri_clean, from = a1, to = a2, 
                                          n = n_sampling, bw = h0, kernel = "gaussian")
                    KDEs[, i] <- kde_result$y
                  }, error = function(e) {
                    KDEs[, i] <- rep(0, n_sampling)
                  })
                } else {
                  KDEs[, i] <- rep(0, n_sampling)
                }
              } else {
                KDEs[, i] <- rep(0, n_sampling)
              }
            }
            
            # Calculate Jensen-Shannon distances - exact translation
            for (t in 1:max_cat_r) {
              for (h in t:max_cat_r) {
                tryCatch({
                  jsd <- gudmm_jensen_shannon(KDEs[, h], KDEs[, t], base = 2)
                  
                  if (is.na(jsd) || jsd > 1) {
                    dist_r_s_th[[r - no_f_cont]][[s]][t, h] <- 0
                  } else {
                    dist_r_s_th[[r - no_f_cont]][[s]][t, h] <- jsd
                  }
                }, error = function(e) {
                  dist_r_s_th[[r - no_f_cont]][[s]][t, h] <- 0
                })
                
                dist_r_s_th[[r - no_f_cont]][[s]][h, t] <- dist_r_s_th[[r - no_f_cont]][[s]][t, h]
              }
            }
          }
        }
      }
      
      # Calculate phi_rth for categorical-categorical interactions - exact translation
      for (t in 1:max_cat_r) {
        for (h in (t+1):max_cat_r) {
          if (h <= max_cat_r) {
            phi_rth_s <- numeric(no_f_disc)
            
            for (s in (no_f_cont + 1):no_f) {
              if (s <= no_f && R[r, s] > 0) {
                s_idx <- s - no_f_cont
                if (s_idx <= no_f_disc && s_idx > 0) {
                  max_cat_s <- max_cats[s_idx]
                  P_rth_ss <- numeric(max_cat_s)
                  E_rth_ss <- numeric(max_cat_s)
                  
                  for (ss in 1:max_cat_s) {
                    # Exact probability calculation as in Python
                    cond1 <- !is.na(X[, r]) & !is.na(X[, s]) & X[, r] == t & X[, s] == ss
                    cond2 <- !is.na(X[, r]) & !is.na(X[, s]) & X[, r] == h & X[, s] == ss
                    
                    P_rth_ss[ss] <- (sum(cond1) + sum(cond2)) / N
                    
                    if (P_rth_ss[ss] > 0) {
                      E_rth_ss[ss] <- -P_rth_ss[ss] * log2(P_rth_ss[ss])
                    }
                  }
                  
                  phi_rth_s[s_idx] <- sum(E_rth_ss) / S_As[s_idx]
                }
              }
            }
            
            # Calculate weighted sum - exact translation
            if (no_f_cont < no_f && length(phi_rth_s) > 0) {
              R_indices <- (no_f_cont + 1):no_f
              R_values <- R[r, R_indices]
              R_E_rth_s <- phi_rth_s * R_values
              phi_rth[[r - no_f_cont]][t, h] <- sum(R_E_rth_s)
              phi_rth[[r - no_f_cont]][h, t] <- phi_rth[[r - no_f_cont]][t, h]
            }
          }
        }
      }
    }
  }
  
  # Combination of distances - exact translation
  for (r in 1:no_f_disc) {
    max_cat_r <- max(X[, r + no_f_cont], na.rm = TRUE)
    
    for (t in 1:max_cat_r) {
      for (h in t:max_cat_r) {
        d_sr_th <- 0
        
        # Sum over continuous features
        for (s in 1:no_f_cont) {
          d_sr_th <- d_sr_th + dist_r_s_th[[r]][[s]][t, h] * R[r + no_f_cont, s]
        }
        
        dist_rth[[r]][t, h] <- d_sr_th
        dist_rth[[r]][h, t] <- d_sr_th
        
        # Final distance calculation - exact as Python
        if (DM == 'DM1' || DM == 'DM2') {
          D_rth[[r]][t, h] <- phi_rth[[r]][t, h] / no_f_disc
        } else {
          D_rth[[r]][t, h] <- (phi_rth[[r]][t, h] + dist_rth[[r]][t, h]) / no_f
        }
        
        D_rth[[r]][h, t] <- D_rth[[r]][t, h]
      }
    }
  }
  
  return(D_rth)
}
