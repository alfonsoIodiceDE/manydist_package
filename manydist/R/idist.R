# requires: purrr (for map)
idist <- function(D, cat_data, pi_nn = 0.1, decision = "prior_corrected") {
  stopifnot(is.matrix(D), nrow(D) == ncol(D))
  n <- nrow(D)
  
  # 1) Precompute neighbor order ONCE for this D
  nn_order <- {
    mat <- matrix(NA_integer_, n, n - 1)
    for (i in seq_len(n)) {
      ord <- order(D[, i], na.last = NA)
      nn  <- ord[ord != i]
      mat[i, ] <- nn
    }
    mat
  }
  # print(nn_order)
  # 2) Per categorical variable: Δ, build Z, form Z Δ Zᵀ, scale by mean off-diagonal
  by_var_idist <- purrr::map(cat_data, ~ {
    lab <- droplevels(factor(.x))
    if (nlevels(lab) < 2) return(matrix(0, n, n))
    
    idelta <- delta_knn_ba(D = D, labels = lab,
                           pi_nn = pi_nn, decision = decision,
                           nn_order = nn_order)
    
# print(idelta)
    
    Z <- model.matrix(~ lab - 1)              # n × G one-hot
    b_v_d <- Z %*% idelta %*% t(Z)
    
    # scale so mean off-diagonal = 1 (with guard)
    den <- mean(as.dist(b_v_d))
    if (!is.finite(den) || den <= 0) den <- 1
    b_v_d / den
  })
  
  # 3) Sum over categorical variables
  distance <- Reduce(`+`, by_var_idist)
  
  # Symmetrize and zero diag for safety
  distance <- (distance + t(distance)) / 2
  diag(distance) <- 0
  distance
}
