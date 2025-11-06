gudmm_mutual_info_score <- function(x, y) {
  tryCatch({
    return(mi.empirical(table(x[, 1], y[, 1])))
  }, error = function(e) return(0))
}
