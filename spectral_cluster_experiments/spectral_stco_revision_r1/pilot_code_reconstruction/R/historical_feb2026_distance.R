# Standalone reconstruction of the mixed-data custom branch at Git commit
# 79563a7 (3 February 2026): Euclidean distance on normalized full PC scores
# plus the unweighted categorical total-variation-profile distance.
#
# For the six numeric columns in Simulation 2, full PCA is an orthogonal
# rotation, so Euclidean distances on PC scores equal Euclidean distances on
# the standardized original columns. The block below evaluates that equivalent
# expression without depending on a historical recipes installation.
#
# The historical mdist() custom Euclidean branch did not consult interaction.
# Its signature also did not accept prop_nn or score. Accordingly this is a
# comparator for the February source, not proof of the exact archived run.

historical_feb2026_sim2_distance <- function(data, components = FALSE) {
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  numeric_data <- data[vapply(data, is.numeric, logical(1))]
  categorical_data <- data[vapply(data, is.factor, logical(1))]
  if (ncol(numeric_data) < 2L || ncol(categorical_data) < 2L) {
    stop("At least two numeric and two factor columns are required.", call. = FALSE)
  }
  if (anyNA(data)) stop("Missing values are not supported.", call. = FALSE)

  standardized <- scale(as.matrix(numeric_data), center = TRUE, scale = TRUE)
  if (any(!is.finite(standardized))) {
    stop("The numeric columns must have positive finite sample SD.", call. = FALSE)
  }
  numeric_distance <- as.matrix(stats::dist(standardized, method = "euclidean"))

  # Recreate historical z_preproc(): one-hot columns in factor/level order,
  # cross-product counts, zero diagonal, then divide each row by its count.
  # R's historical `ZZ / zm` recycles `zm` down columns, hence it scales rows.
  categorical_data <- lapply(categorical_data, droplevels)
  level_counts <- vapply(categorical_data, nlevels, integer(1))
  indicator <- do.call(cbind, lapply(categorical_data, function(column) {
    vapply(levels(column), function(level) as.numeric(column == level),
           numeric(length(column)))
  }))
  count <- colSums(indicator)
  if (any(count == 0L)) stop("An empty factor level remains.", call. = FALSE)
  profile <- crossprod(indicator)
  diag(profile) <- 0
  profile <- sweep(profile, 1L, count, FUN = "/")

  # historical cat_custom_delta("tot_var_dist") takes the Manhattan distance
  # between full category profiles, divides by 2 * (q - 1), and zeros every
  # off-diagonal variable block. cdist() then computes Z delta Z'.
  q <- length(level_counts)
  delta <- matrix(0, ncol(indicator), ncol(indicator))
  stop_at <- cumsum(level_counts)
  start_at <- c(1L, head(stop_at, -1L) + 1L)
  for (variable in seq_len(q)) {
    index <- start_at[[variable]]:stop_at[[variable]]
    delta[index, index] <- as.matrix(stats::dist(
      profile[index, , drop = FALSE], method = "manhattan"
    )) / (2 * (q - 1L))
  }
  categorical_distance <- indicator %*% delta %*% t(indicator)
  distance <- numeric_distance + categorical_distance
  if (!isTRUE(components)) return(distance)
  list(
    distance = distance,
    numeric = numeric_distance,
    categorical = categorical_distance,
    delta = delta,
    level_counts = level_counts
  )
}
