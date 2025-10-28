gudmm_mutual_info_regression <- function(x, y) {
  tryCatch({
    x_disc <- cut(x[, 1], breaks = min(10, length(unique(x[, 1]))), labels = FALSE)
    y_disc <- cut(y[, 1], breaks = min(10, length(unique(y[, 1]))), labels = FALSE)
    if (any(is.na(x_disc)) || any(is.na(y_disc))) return(0)
    return(mi.empirical(table(x_disc, y_disc)))
  }, error = function(e) return(0))
}
