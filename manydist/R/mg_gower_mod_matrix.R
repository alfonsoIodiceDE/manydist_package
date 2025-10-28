mg_gower_mod_matrix <- function(x, use_weights = TRUE) {
  source("R/mg_gower_dist_modify.R")
  source("R/mg_self_adaptive_distance.R")
  source("R/mg_normalized_MI.R")
  
  # Convert tibble to data.frame to ensure proper indexing behavior
  if (is_tibble(x)) {
    x <- as.data.frame(x)
  }
  
  num_idx <- sapply(x, is.numeric)
  X_num   <- x[,  num_idx, drop = FALSE]
  X_cat   <- x[, !num_idx, drop = FALSE]
  
  if (use_weights) {
    # Get feature importance
    feature_importance1 <- mg_self_adaptive_distance(X_num, X_cat)
    feature_importance2 <- mg_normalized_MI(X_num, X_cat)
    feature_importance <- feature_importance1 * feature_importance2
    feature_importance <- feature_importance / sum(feature_importance)
  } else {
    feature_importance <- rep(1, ncol(x))
  }
  
  mg_gower_dist_modify(data.x = x,
                    var.weights = feature_importance,
                    robcb = "iqr",          
                    KR.corr = TRUE)
}
