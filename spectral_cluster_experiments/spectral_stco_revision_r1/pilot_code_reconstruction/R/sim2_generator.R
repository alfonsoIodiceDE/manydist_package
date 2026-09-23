# Reconstruction of the available historical Simulation 2 script.
# This is not asserted to be the exact driver of the archived 50-replicate data.

sim2_reconstructed_categories <- function(n, index_rule = c("literal", "corrected")) {
  index_rule <- match.arg(index_rule)
  if (length(n) != 1L || !is.finite(n) || n != as.integer(n) ||
      n < 100L || n %% 10L != 0L) {
    stop("n must be an integer >= 100 divisible by 10.", call. = FALSE)
  }
  n <- as.integer(n)

  # srswor() in the historical script selects a uniformly sampled fixed-size
  # subset. Base R gives the same sampling design, but not the same RNG stream.
  chosen_1 <- sort(sample.int(n, size = as.integer(0.6 * n)))
  category_1 <- rep(1L, n)
  if (identical(index_rule, "literal")) {
    # Preserve the historical script's unparenthesized R index expressions.
    category_1[chosen_1[1:n * 0.2]] <- 2L
    category_1[chosen_1[(n * 0.2 + 1):n * 0.6]] <- 3L
  } else {
    category_1[chosen_1[seq_len(as.integer(0.2 * n))]] <- 2L
    category_1[chosen_1[seq.int(as.integer(0.2 * n + 1), as.integer(0.6 * n))]] <- 3L
  }

  category_2 <- rep(1L, n)
  category_2[sample.int(n, size = as.integer(0.2 * n))] <- 2L

  chosen_4 <- sort(sample.int(n, size = as.integer(0.6 * n)))
  category_4 <- rep(1L, n)
  if (identical(index_rule, "literal")) {
    category_4[chosen_4[1:n * 0.3]] <- 2L
    category_4[chosen_4[(n * 0.3 + 1):n * 0.6]] <- 3L
  } else {
    category_4[chosen_4[seq_len(as.integer(0.3 * n))]] <- 2L
    category_4[chosen_4[seq.int(as.integer(0.3 * n + 1), as.integer(0.6 * n))]] <- 3L
  }

  category_3 <- rep(1L, n)
  category_3[sample.int(n, size = as.integer(0.4 * n))] <- 2L

  list(
    C1 = factor(category_1, levels = 1:3),
    C2 = factor(category_2, levels = 1:2),
    C3 = factor(category_3, levels = 1:2),
    C4 = factor(category_4, levels = 1:3)
  )
}

sim2_reconstructed_data <- function(
    seed,
    n = 500L,
    index_rule = c("literal", "corrected")
) {
  index_rule <- match.arg(index_rule)
  if (length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be a finite integer.", call. = FALSE)
  }
  set.seed(as.integer(seed))

  # The historical script used SimDesign::rmvnorm with identity covariance.
  # Base R reproduces that distribution, not the package-specific RNG stream.
  numeric_data <- matrix(stats::rnorm(n * 6L), nrow = n, ncol = 6L)
  categories <- sim2_reconstructed_categories(n, index_rule)
  coefficient_matrix <- matrix(
    c(0.8, 0.5, -0.4, 0.8, -0.4, 0.5),
    nrow = 3L,
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(levels(categories$C1), levels(categories$C2))
  )
  coefficient <- coefficient_matrix[cbind(
    as.integer(categories$C1),
    as.integer(categories$C2)
  )]
  intercept <- ifelse(coefficient == 0.5, 8, 0)

  for (variable in 3:1) {
    numeric_data[, variable] <- intercept + coefficient *
      rowSums(numeric_data[, (variable + 1L):6L, drop = FALSE])
  }
  colnames(numeric_data) <- paste0("V", seq_len(6L))
  truth <- factor(
    ifelse(coefficient == 0.8, "C1", ifelse(coefficient == -0.4, "C2", "C3")),
    levels = c("C1", "C2", "C3")
  )
  data <- data.frame(numeric_data, categories, check.names = FALSE)

  list(
    x = data,
    truth = truth,
    diagnostics = list(
      source = "available historical simulation xBIG.R, not verified archived driver",
      index_rule = index_rule,
      n = n,
      coefficient_matrix = coefficient_matrix,
      intercept_for_coefficient_0.5 = 8,
      category_counts = lapply(categories, table),
      truth_counts = table(truth),
      rng_caveat = paste(
        "Base-R normal and fixed-subset draws match the historical statistical",
        "design but not SimDesign/sampling seed-for-seed streams."
      )
    )
  )
}
