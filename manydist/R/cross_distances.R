# Coordinate differences avoid both backend-specific failures and cancellation
# from subtracting nearly equal squared norms. Rows are test-to-training.
.cross_euclidean_distance <- function(xnew, x) {
  xnew <- as.matrix(xnew)
  x <- as.matrix(x)
  if (ncol(xnew) != ncol(x)) {
    stop("Training and new data must have the same number of coordinates.",
         call. = FALSE)
  }
  squared <- matrix(0, nrow = nrow(xnew), ncol = nrow(x))
  for (column in seq_len(ncol(x))) {
    squared <- squared + outer(xnew[, column], x[, column], "-")^2
  }
  sqrt(squared)
}
