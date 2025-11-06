gudmm_mutual_info_classif <- function(x, y) {
  tryCatch({
    x_disc <- cut(x[, 1], breaks = min(10, length(unique(x[, 1]))), labels = FALSE)
    if (any(is.na(x_disc))) return(0)
    return(mi.empirical(table(x_disc, y[, 1])))
  }, error = function(e) return(0))
}
