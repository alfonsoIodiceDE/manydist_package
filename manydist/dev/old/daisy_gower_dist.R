daisy_gower_dist <- function(train, new_data, ...) {
  X <- dplyr::select(train, -Species)
  Z <- dplyr::select(new_data, -Species)

  full   <- rbind(Z, X)
  Dfull  <- as.matrix(cluster::daisy(full, metric = "gower"))

  n_test  <- nrow(Z)
  n_train <- nrow(X)

  # rectangular test × train block
  Dfull[1:n_test, (n_test + 1):(n_test + n_train), drop = FALSE]
}
