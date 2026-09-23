full_assert_scalar <- function(x, label, lower = -Inf, upper = Inf, integer = FALSE) {
  if (length(x) != 1L || !is.finite(x) || x < lower || x > upper ||
      (isTRUE(integer) && x != as.integer(x))) {
    stop("Invalid `", label, "` in scenario configuration.", call. = FALSE)
  }
  invisible(x)
}

full_ar1_covariance <- function(p, rho) {
  full_assert_scalar(p, "p", lower = 1, integer = TRUE)
  full_assert_scalar(rho, "ar1_rho", lower = -0.99, upper = 0.99)
  sigma <- rho^abs(outer(seq_len(p), seq_len(p), `-`))
  min_eigenvalue <- min(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(min_eigenvalue) || min_eigenvalue <= 1e-10) {
    stop("The configured covariance matrix is not positive definite.", call. = FALSE)
  }
  # Cholesky is an additional executable guard against invalid covariance.
  chol(sigma)
  list(matrix = sigma, min_eigenvalue = min_eigenvalue)
}

full_balanced_labels <- function(n, levels) {
  factor(
    sample(rep(levels, length.out = n), size = n, replace = FALSE),
    levels = levels
  )
}

full_marginal_factor <- function(truth, variable_id, accuracy = 0.80) {
  truth <- as.integer(factor(truth))
  k <- max(truth)
  levels <- paste0("L", seq_len(k))
  preferred <- ((truth + variable_id - 2L) %% k) + 1L
  observed <- preferred
  mistakes <- stats::runif(length(truth)) > accuracy
  for (i in which(mistakes)) {
    observed[i] <- sample(setdiff(seq_len(k), preferred[i]), 1L)
  }
  factor(levels[observed], levels = levels)
}

full_contaminate_factor <- function(x, rate) {
  x <- droplevels(factor(x))
  full_assert_scalar(rate, "categorical_contamination", lower = 0, upper = 1)
  n_change <- as.integer(round(rate * length(x)))
  if (!n_change) {
    return(list(value = x, changed = 0L, realized_rate = 0))
  }
  positions <- sample(seq_along(x), n_change, replace = FALSE)
  integer_values <- as.integer(x)
  for (position in positions) {
    integer_values[position] <- sample(
      setdiff(seq_len(nlevels(x)), integer_values[position]),
      1L
    )
  }
  list(
    value = factor(levels(x)[integer_values], levels = levels(x)),
    changed = n_change,
    realized_rate = n_change / length(x)
  )
}

full_innovations <- function(n, p, family) {
  if (identical(family, "gaussian")) {
    return(matrix(stats::rnorm(n * p), nrow = n, ncol = p))
  }
  if (family %in% c("student_t3", "student_t5")) {
    df <- if (identical(family, "student_t3")) 3 else 5
    gaussian <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    radial <- sqrt(stats::rchisq(n, df = df) / df)
    # The multiplier makes each marginal variance one before correlation.
    return((gaussian / radial) * sqrt((df - 2) / df))
  }
  if (identical(family, "centered_chisq3")) {
    return(matrix(
      (stats::rchisq(n * p, df = 3) - 3) / sqrt(6),
      nrow = n,
      ncol = p
    ))
  }
  stop("Unknown continuous family: ", family, call. = FALSE)
}

full_make_controlled_overlap <- function(parameters, seed) {
  n <- as.integer(parameters$n)
  k <- as.integer(parameters$k)
  p_signal <- as.integer(parameters$p_numeric_signal)
  p_noise <- as.integer(parameters$p_numeric_noise %||% 0L)
  q_signal <- as.integer(parameters$q_categorical_signal)
  q_noise <- as.integer(parameters$q_categorical_noise %||% 0L)
  categorical_accuracy <- as.numeric(parameters$categorical_accuracy)
  target_overlap <- as.numeric(parameters$adjacent_overlap)
  family <- as.character(parameters$family)
  rho <- as.numeric(parameters$ar1_rho)

  full_assert_scalar(n, "n", lower = 30, integer = TRUE)
  full_assert_scalar(k, "k", lower = 2, integer = TRUE)
  full_assert_scalar(p_signal, "p_numeric_signal", lower = 2, integer = TRUE)
  full_assert_scalar(p_noise, "p_numeric_noise", lower = 0, integer = TRUE)
  full_assert_scalar(q_signal, "q_categorical_signal", lower = 1, integer = TRUE)
  full_assert_scalar(q_noise, "q_categorical_noise", lower = 0, integer = TRUE)
  full_assert_scalar(categorical_accuracy, "categorical_accuracy", lower = 0.34, upper = 1)
  full_assert_scalar(target_overlap, "adjacent_overlap", lower = 0.01, upper = 0.80)

  covariance <- full_ar1_covariance(p_signal, rho)
  sigma <- covariance$matrix
  sigma_inverse <- solve(sigma)
  mahalanobis_step <- -2 * stats::qnorm(target_overlap / 2)
  raw_step <- mahalanobis_step / sqrt(sigma_inverse[1L, 1L])
  centers <- ((seq_len(k) - (k + 1) / 2) * raw_step)

  # Keep independent, deterministic streams for the truth, continuous signal,
  # numerical noise, categorical signal, and categorical noise.  This makes
  # the irrelevant-variable grids genuinely paired: adding noise variables
  # cannot silently change any of the fixed signal variables.
  set.seed(seed)
  truth_integer <- sample(rep(seq_len(k), length.out = n), n, replace = FALSE)
  signal <- matrix(NA_real_, nrow = n, ncol = p_signal)
  chol_sigma <- chol(sigma)
  set.seed(seed + 1L)
  for (cluster in seq_len(k)) {
    rows <- which(truth_integer == cluster)
    innovation <- full_innovations(length(rows), p_signal, family)
    shifted <- innovation %*% chol_sigma
    shifted[, 1L] <- shifted[, 1L] + centers[cluster]
    signal[rows, ] <- shifted
  }
  colnames(signal) <- paste0("V_signal_", seq_len(p_signal))

  numeric_noise <- if (p_noise) {
    set.seed(seed + 2L)
    value <- matrix(stats::rnorm(n * p_noise), nrow = n, ncol = p_noise)
    colnames(value) <- paste0("V_noise_", seq_len(p_noise))
    value
  } else {
    NULL
  }

  truth <- factor(truth_integer, levels = seq_len(k), labels = paste0("C", seq_len(k)))
  set.seed(seed + 3L)
  categorical <- lapply(seq_len(q_signal), function(j) {
    full_marginal_factor(truth, j, categorical_accuracy)
  })
  names(categorical) <- paste0("C_signal_", seq_len(q_signal))
  if (q_noise) {
    set.seed(seed + 4L)
    noise_factors <- lapply(seq_len(q_noise), function(j) {
      full_balanced_labels(n, c("N1", "N2", "N3"))
    })
    names(noise_factors) <- paste0("C_noise_", seq_len(q_noise))
    categorical <- c(categorical, noise_factors)
  }

  x <- data.frame(signal, check.names = FALSE)
  if (!is.null(numeric_noise)) {
    x <- data.frame(x, numeric_noise, check.names = FALSE)
  }
  x <- data.frame(x, as.data.frame(categorical), check.names = FALSE)
  list(
    x = x,
    truth = truth,
    diagnostics = list(
      generator = "controlled_overlap",
      family = family,
      n = n,
      k = k,
      p_numeric_signal = p_signal,
      p_numeric_noise = p_noise,
      q_categorical_signal = q_signal,
      q_categorical_noise = q_noise,
      ar1_rho = rho,
      covariance_min_eigenvalue = covariance$min_eigenvalue,
      target_adjacent_overlap = target_overlap,
      adjacent_mahalanobis_distance = mahalanobis_step,
      overlap_interpretation = if (family == "gaussian") {
        "exact two-component equal-covariance overlap coefficient for adjacent centers"
      } else {
        "Gaussian calibration retained for geometry; not an exact non-Gaussian overlap coefficient"
      }
    )
  )
}

full_assign_balanced_within <- function(groups, target_levels) {
  groups <- factor(groups)
  out <- character(length(groups))
  for (group in levels(groups)) {
    rows <- which(groups == group)
    out[rows] <- sample(
      rep(target_levels, length.out = length(rows)),
      replace = FALSE
    )
  }
  factor(out, levels = target_levels)
}

full_make_moons <- function(parameters, seed) {
  set.seed(seed)
  n_per_moon <- as.integer(parameters$n_per_moon)
  state <- as.character(parameters$interaction_state)
  contamination_rate <- as.numeric(parameters$categorical_contamination %||% 0)
  q_noise <- as.integer(parameters$q_categorical_noise %||% 0L)
  noise_sd <- as.numeric(parameters$moon_noise_sd)

  full_assert_scalar(n_per_moon, "n_per_moon", lower = 18, integer = TRUE)
  if (n_per_moon %% 9L != 0L) {
    stop("`n_per_moon` must be divisible by 9 for exact balancing.", call. = FALSE)
  }
  if (!state %in% c("signal", "null")) {
    stop("`interaction_state` must be `signal` or `null`.", call. = FALSE)
  }

  n_per_band <- n_per_moon / 3L
  draw_stratified_theta <- function() {
    unlist(lapply(0:2, function(band) {
      stats::runif(
        n_per_band,
        min = band * pi / 3,
        max = (band + 1) * pi / 3
      )
    }), use.names = FALSE)
  }
  theta_upper <- draw_stratified_theta()
  theta_lower <- draw_stratified_theta()
  upper <- data.frame(
    V1 = cos(theta_upper) + stats::rnorm(n_per_moon, sd = noise_sd),
    V2 = sin(theta_upper) + stats::rnorm(n_per_moon, sd = noise_sd),
    arc = theta_upper / pi,
    moon = "upper"
  )
  lower <- data.frame(
    V1 = 1 - cos(theta_lower) + stats::rnorm(n_per_moon, sd = noise_sd),
    V2 = 0.50 - sin(theta_lower) + stats::rnorm(n_per_moon, sd = noise_sd),
    arc = theta_lower / pi,
    moon = "lower"
  )
  data <- rbind(upper, lower)
  data$moon <- factor(data$moon, levels = c("upper", "lower"))
  equal_thirds <- rep(c("A", "B", "C"), each = n_per_moon / 3L)

  clean_band <- rep(equal_thirds, times = 2L)
  if (state == "null") {
    for (moon in levels(data$moon)) {
      rows <- which(data$moon == moon)
      clean_band[rows] <- sample(equal_thirds)
    }
  }
  clean_band <- factor(clean_band, levels = c("A", "B", "C"))

  # Randomize row order after the stratified draws so every replicate is a
  # genuinely resampled dataset while retaining exact band counts.
  permutation <- sample(seq_len(nrow(data)))
  data <- data[permutation, , drop = FALSE]
  clean_band <- clean_band[permutation]

  latent_truth <- if (state == "signal") {
    interaction(data$moon, clean_band, sep = " x ", drop = TRUE)
  } else {
    data$moon
  }
  nuisance <- full_assign_balanced_within(
    interaction(data$moon, clean_band, drop = TRUE),
    c("N1", "N2", "N3")
  )
  contaminated <- full_contaminate_factor(clean_band, contamination_rate)

  categorical <- list(band = contaminated$value, nuisance = nuisance)
  if (q_noise) {
    extra_noise <- lapply(seq_len(q_noise), function(j) {
      full_assign_balanced_within(latent_truth, c("N1", "N2", "N3"))
    })
    names(extra_noise) <- paste0("cat_noise_", seq_len(q_noise))
    categorical <- c(categorical, extra_noise)
  }

  x <- data.frame(data[c("V1", "V2")], categorical, check.names = FALSE)
  list(
    x = x,
    truth = droplevels(factor(latent_truth)),
    diagnostics = list(
      generator = "moons",
      interaction_state = state,
      n = nrow(x),
      n_per_moon = n_per_moon,
      moon_noise_sd = noise_sd,
      requested_contamination = contamination_rate,
      contaminated_count = contaminated$changed,
      realized_contamination = contaminated$realized_rate,
      q_categorical_noise = q_noise
    )
  )
}

full_submitted_equicorrelation <- function(p, rho) {
  full_assert_scalar(p, "p_numeric_signal", lower = 1, integer = TRUE)
  full_assert_scalar(rho, "rho", lower = -1, upper = 1)

  requested <- matrix(rho, nrow = p, ncol = p)
  diag(requested) <- 1
  decomposition <- eigen(requested, symmetric = TRUE)
  retained_values <- pmax(decomposition$values, 0)
  root <- decomposition$vectors %*%
    diag(sqrt(retained_values), nrow = p) %*%
    t(decomposition$vectors)

  list(
    root = root,
    requested_min_eigenvalue = min(decomposition$values),
    clipped_eigenvalues = sum(decomposition$values < 0),
    realized_rank = sum(retained_values > sqrt(.Machine$double.eps))
  )
}

full_make_submitted_simulation_1 <- function(parameters, seed) {
  n <- as.integer(parameters$n)
  p_signal <- as.integer(parameters$p_numeric_signal)
  p_noise <- as.integer(parameters$p_numeric_noise %||% 0L)
  q <- as.integer(parameters$q_categorical)

  full_assert_scalar(n, "n", lower = 30, integer = TRUE)
  full_assert_scalar(p_signal, "p_numeric_signal", lower = 1, integer = TRUE)
  full_assert_scalar(p_noise, "p_numeric_noise", lower = 0, integer = TRUE)
  full_assert_scalar(q, "q_categorical", lower = 1, integer = TRUE)

  cluster_sizes <- c(round(0.20 * n), round(0.30 * n))
  cluster_sizes <- c(cluster_sizes, n - sum(cluster_sizes))
  cluster_means <- list(
    rep(0, p_signal),
    rep(4, p_signal),
    rep(c(-1, 5), length.out = p_signal)
  )
  cluster_correlations <- c(0.50, 0.80, -0.50)
  categorical_slopes <- c(2, -2, 2)

  set.seed(seed)
  numeric_signal <- matrix(NA_real_, nrow = n, ncol = p_signal)
  categorical <- matrix(NA_integer_, nrow = n, ncol = q)
  truth_integer <- rep(seq_len(3L), times = cluster_sizes)
  covariance_diagnostics <- vector("list", 3L)

  row_start <- 1L
  for (cluster in seq_len(3L)) {
    rows <- seq.int(row_start, length.out = cluster_sizes[[cluster]])
    row_start <- row_start + cluster_sizes[[cluster]]
    covariance <- full_submitted_equicorrelation(
      p_signal,
      cluster_correlations[[cluster]]
    )
    block <- matrix(
      stats::rnorm(length(rows) * p_signal),
      nrow = length(rows),
      ncol = p_signal
    ) %*% covariance$root
    block <- sweep(block, 2L, cluster_means[[cluster]], `+`)
    numeric_signal[rows, ] <- block

    for (variable in seq_len(q)) {
      previous_sum <- if (variable == 1L) {
        rep(0, length(rows))
      } else {
        rowSums(categorical[rows, seq_len(variable - 1L), drop = FALSE])
      }
      linear_predictor <- 1 + categorical_slopes[[cluster]] *
        (rowSums(block) + previous_sum)
      probability_one <- stats::plogis(-linear_predictor)
      categorical[rows, variable] <- stats::rbinom(
        length(rows),
        size = 1L,
        prob = probability_one
      )
    }
    covariance_diagnostics[[cluster]] <- list(
      cluster = cluster,
      requested_rho = cluster_correlations[[cluster]],
      requested_min_eigenvalue = covariance$requested_min_eigenvalue,
      clipped_eigenvalues = covariance$clipped_eigenvalues,
      realized_rank = covariance$realized_rank
    )
  }

  colnames(numeric_signal) <- paste0("V", seq_len(p_signal))
  numeric_noise <- if (p_noise) {
    value <- matrix(stats::runif(n * p_noise), nrow = n, ncol = p_noise)
    colnames(value) <- paste0("V_noise_", seq_len(p_noise))
    value
  } else {
    NULL
  }
  categorical <- as.data.frame(lapply(seq_len(q), function(variable) {
    factor(categorical[, variable], levels = c(0, 1))
  }))
  names(categorical) <- paste0("C", seq_len(q))

  x <- data.frame(numeric_signal, check.names = FALSE)
  if (!is.null(numeric_noise)) {
    x <- data.frame(x, numeric_noise, check.names = FALSE)
  }
  x <- data.frame(x, categorical, check.names = FALSE)

  list(
    x = x,
    truth = factor(truth_integer, levels = seq_len(3L), labels = paste0("C", seq_len(3L))),
    diagnostics = list(
      generator = "submitted_simulation_1",
      source = "submitted_manuscript",
      n = n,
      cluster_proportions = c(0.20, 0.30, 0.50),
      p_numeric_signal = p_signal,
      p_numeric_noise = p_noise,
      q_categorical = q,
      covariance = covariance_diagnostics,
      covariance_caveat = paste(
        "The requested rho = -0.5 equicorrelation matrix is indefinite",
        "when p > 3; negative eigenvalues are clipped to zero to reproduce",
        "the effective eigen-based legacy draw without silently replacing",
        "the submitted specification."
      ),
      historical_code_audit = paste(
        "This generator follows the submitted manuscript. The historical",
        "script constructed the categorical variables from unshifted latent",
        "numerics, used different predictor indexing, and used n = 500/1000",
        "rather than the manuscript's n = 1000/2000."
      )
    )
  )
}

full_fixed_proportion_factor <- function(n, probabilities, labels) {
  if (length(probabilities) != length(labels) ||
      any(!is.finite(probabilities)) || any(probabilities < 0) ||
      abs(sum(probabilities) - 1) > sqrt(.Machine$double.eps)) {
    stop("Invalid fixed factor proportions.", call. = FALSE)
  }
  counts <- as.integer(round(n * probabilities))
  counts[[length(counts)]] <- n - sum(counts[-length(counts)])
  factor(
    sample(rep(labels, times = counts), size = n, replace = FALSE),
    levels = labels
  )
}

full_make_submitted_simulation_2 <- function(parameters, seed) {
  n <- as.integer(parameters$n)
  full_assert_scalar(n, "n", lower = 30, integer = TRUE)

  set.seed(seed)
  c1 <- full_fixed_proportion_factor(n, c(0.40, 0.20, 0.40), c("1", "2", "3"))
  c2 <- full_fixed_proportion_factor(n, c(0.80, 0.20), c("1", "2"))
  c3 <- full_fixed_proportion_factor(n, c(0.60, 0.40), c("1", "2"))

  coefficient_matrix <- matrix(
    c(0.8, -0.4, -0.4, 0.1, 0.8, 0.1),
    nrow = 3L,
    ncol = 2L,
    dimnames = list(levels(c1), levels(c2))
  )
  coefficient <- coefficient_matrix[cbind(as.integer(c1), as.integer(c2))]
  intercept <- ifelse(coefficient == 0.1, 3, 0)

  numeric_data <- matrix(NA_real_, nrow = n, ncol = 6L)
  numeric_data[, 4:6] <- matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L)
  for (variable in 3:1) {
    numeric_data[, variable] <- intercept + coefficient *
      rowSums(numeric_data[, (variable + 1L):6L, drop = FALSE])
  }
  colnames(numeric_data) <- paste0("V", seq_len(6L))

  truth <- factor(
    ifelse(coefficient == 0.8, "C1", ifelse(coefficient == -0.4, "C2", "C3")),
    levels = c("C1", "C2", "C3")
  )
  x <- data.frame(numeric_data, C1 = c1, C2 = c2, C3 = c3, check.names = FALSE)

  list(
    x = x,
    truth = truth,
    diagnostics = list(
      generator = "submitted_simulation_2",
      source = "submitted_manuscript",
      n = n,
      p_numeric = 6L,
      q_categorical = 3L,
      coefficient_matrix = coefficient_matrix,
      intercept_for_coefficient_0.1 = 3,
      historical_code_audit = paste(
        "This generator follows the submitted manuscript. The separate",
        "historical scripts used 0.5 instead of 0.1, an intercept of 8",
        "instead of 3, and four rather than three categorical variables."
      )
    )
  )
}

# Code-faithful reconstructions of the data generators used by the historical
# Simulation 1 scripts.  These deliberately retain two non-obvious features of
# those scripts: categorical variables are generated from the unshifted
# numerical draws, and categorical variable j uses the historical
# 1:(j + p_numeric - 1) predictor slice (so C1 omits the final numerical
# variable).  The scripts cbind the initially zero categorical matrix into
# `XX`; later writes to `Xcat` do not update that copy.  Consequently, the
# apparent lagged categorical predictors remain zero and C2,...,Cq use the
# same numerical linear predictor.  Keep these generators separate from the
# manuscript-specification generators above: they answer different audit
# questions and neither should silently replace the other.
full_make_legacy_simulation_1 <- function(parameters, seed) {
  n <- as.integer(parameters$n)
  p_numeric <- as.integer(parameters$p_numeric %||% parameters$p_numeric_signal)
  q_categorical <- as.integer(parameters$q_categorical)

  full_assert_scalar(n, "n", lower = 10, integer = TRUE)
  full_assert_scalar(p_numeric, "p_numeric", lower = 2, integer = TRUE)
  full_assert_scalar(q_categorical, "q_categorical", lower = 1, integer = TRUE)
  required_cluster_sizes <- n * c(0.20, 0.30, 0.50)
  if (any(required_cluster_sizes != as.integer(required_cluster_sizes))) {
    stop(
      "`n` must make 20%, 30%, and 50% integer cluster sizes for the ",
      "historical Simulation 1 construction.",
      call. = FALSE
    )
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("The legacy Simulation 1 generator requires the `mvtnorm` package.",
         call. = FALSE)
  }

  cluster_sizes <- as.integer(required_cluster_sizes)
  correlations <- c(0.50, 0.80, -0.50)
  categorical_slopes <- c(2, -2, 2)
  numeric_shifts <- list(
    rep(0, p_numeric),
    rep(4, p_numeric),
    rep(c(-1, 5), length.out = p_numeric)
  )

  set.seed(seed)
  blocks <- vector("list", 3L)
  covariance_diagnostics <- vector("list", 3L)
  total_variables <- p_numeric + q_categorical

  for (cluster in seq_len(3L)) {
    covariance <- matrix(
      correlations[[cluster]],
      nrow = p_numeric,
      ncol = p_numeric
    )
    diag(covariance) <- 1
    decomposition <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values

    numeric_unshifted <- suppressWarnings(mvtnorm::rmvnorm(
      cluster_sizes[[cluster]],
      mean = rep(0, p_numeric),
      sigma = covariance,
      method = "eigen"
    ))

    categorical <- matrix(
      0,
      nrow = cluster_sizes[[cluster]],
      ncol = q_categorical
    )
    predictor_bank <- cbind(1, numeric_unshifted, categorical)
    beta <- matrix(
      c(1, rep(categorical_slopes[[cluster]], total_variables)),
      nrow = q_categorical,
      ncol = total_variables + 1L,
      byrow = TRUE
    )
    for (row in seq_len(cluster_sizes[[cluster]])) {
      for (variable in seq_len(q_categorical)) {
        columns <- seq_len(variable + p_numeric - 1L)
        linear_predictor <- beta[variable, columns] %*%
          predictor_bank[row, columns]
        categorical[row, variable] <- stats::rbinom(
          1L,
          size = 1L,
          prob = 1 / (1 + exp(linear_predictor))
        )
        # Do not update `predictor_bank`: the legacy `XX <- cbind(XX, Xcat)`
        # copied the all-zero categorical matrix before this loop.
      }
    }

    numeric_shifted <- sweep(
      numeric_unshifted,
      2L,
      numeric_shifts[[cluster]],
      `+`
    )
    blocks[[cluster]] <- cbind(numeric_shifted, categorical)
    covariance_diagnostics[[cluster]] <- list(
      cluster = cluster,
      requested_rho = correlations[[cluster]],
      requested_min_eigenvalue = min(decomposition),
      clipped_eigenvalues = sum(decomposition < 0)
    )
  }

  combined <- do.call(rbind, blocks)
  numeric_names <- paste0("V", seq_len(p_numeric))
  categorical_names <- paste0("C", seq_len(q_categorical))
  x <- data.frame(
    combined[, seq_len(p_numeric), drop = FALSE],
    check.names = FALSE
  )
  names(x) <- numeric_names
  categorical_columns <- lapply(seq_len(q_categorical), function(variable) {
    factor(combined[, p_numeric + variable], levels = c(0, 1))
  })
  names(categorical_columns) <- categorical_names
  x <- data.frame(x, categorical_columns, check.names = FALSE)

  list(
    x = x,
    truth = factor(
      rep(seq_len(3L), times = cluster_sizes),
      levels = seq_len(3L),
      labels = paste0("C", seq_len(3L))
    ),
    diagnostics = list(
      generator = "legacy_simulation_1",
      source = "historical_executable_code",
      source_scripts = c(
        "SetupSim.R",
        "SetupSim_1020.R",
        "SetupSim_2010.R",
        "SetupSim_bignR.R"
      ),
      n = n,
      cluster_proportions = cluster_sizes / n,
      p_numeric = p_numeric,
      q_categorical = q_categorical,
      covariance = covariance_diagnostics,
      categorical_predictor_rule = paste(
        "historical 1:(j + p_numeric - 1) slice with cbind-ed categorical",
        "predictors remaining zero"
      ),
      categorical_generated_before_numeric_shifts = TRUE,
      exact_historical_observations_recoverable = FALSE
    )
  )
}

full_legacy_srswor <- function(sample_size, population_size) {
  sample_size <- as.integer(sample_size)
  population_size <- as.integer(population_size)
  selected <- integer(population_size)
  selected[sample.int(population_size, sample_size)] <- 1L
  selected
}

# Code-faithful reconstruction of `simulation xBIG.R`.  In particular it
# preserves the fourth nuisance categorical variable, beta = 0.5 cell, offset
# 8, sorted-index assignment induced by which(srswor(...)), and recursive
# construction of V3,V2,V1 from the later numerical coordinates.
full_make_legacy_simulation_2 <- function(parameters, seed) {
  n <- as.integer(parameters$n)
  p_numeric <- as.integer(parameters$p_numeric %||% 6L)
  q_categorical <- as.integer(parameters$q_categorical %||% 4L)

  full_assert_scalar(n, "n", lower = 10, integer = TRUE)
  if (p_numeric != 6L || q_categorical != 4L) {
    stop(
      "The historical Simulation 2 code fixes `p_numeric = 6` and ",
      "`q_categorical = 4`.",
      call. = FALSE
    )
  }
  required_fractions <- n * c(0.20, 0.30, 0.40, 0.60)
  if (any(required_fractions != as.integer(required_fractions))) {
    stop(
      "`n` must make 20%, 30%, 40%, and 60% integer counts for the ",
      "historical Simulation 2 construction.",
      call. = FALSE
    )
  }

  set.seed(seed)
  numeric_data <- matrix(stats::rnorm(n * p_numeric), nrow = n)

  c1 <- rep(1L, n)
  selected <- which(full_legacy_srswor(0.60 * n, n) == 1L)
  c1[selected[seq_len(0.20 * n)]] <- 2L
  c1[selected[seq.int(0.20 * n + 1L, 0.60 * n)]] <- 3L

  c2 <- rep(1L, n)
  c2[which(full_legacy_srswor(0.20 * n, n) == 1L)] <- 2L

  c4 <- rep(1L, n)
  selected <- which(full_legacy_srswor(0.60 * n, n) == 1L)
  c4[selected[seq_len(0.30 * n)]] <- 2L
  c4[selected[seq.int(0.30 * n + 1L, 0.60 * n)]] <- 3L

  c3 <- rep(1L, n)
  c3[which(full_legacy_srswor(0.40 * n, n) == 1L)] <- 2L

  coefficient_matrix <- matrix(
    c(0.8, 0.5, -0.4, 0.8, -0.4, 0.5),
    nrow = 3L,
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), c("1", "2"))
  )
  coefficient <- coefficient_matrix[cbind(c1, c2)]
  truth_integer <- rep(3L, n)

  for (variable in 3:1) {
    numeric_data[, variable] <- coefficient * rowSums(
      numeric_data[, (variable + 1L):p_numeric, drop = FALSE]
    )
    numeric_data[coefficient == 0.5, variable] <-
      numeric_data[coefficient == 0.5, variable] + 8
  }
  truth_integer[coefficient == 0.8] <- 1L
  truth_integer[coefficient == -0.4] <- 2L
  colnames(numeric_data) <- paste0("V", seq_len(p_numeric))

  x <- data.frame(
    numeric_data,
    C1 = factor(c1, levels = 1:3),
    C2 = factor(c2, levels = 1:2),
    C3 = factor(c3, levels = 1:2),
    C4 = factor(c4, levels = 1:3),
    check.names = FALSE
  )

  list(
    x = x,
    truth = factor(
      truth_integer,
      levels = seq_len(3L),
      labels = paste0("C", seq_len(3L))
    ),
    diagnostics = list(
      generator = "legacy_simulation_2",
      source = "historical_executable_code",
      source_script = "simulation xBIG.R",
      n = n,
      p_numeric = p_numeric,
      q_categorical = q_categorical,
      coefficient_matrix = coefficient_matrix,
      intercept_for_coefficient_0.5 = 8,
      nuisance_categorical = c("C3", "C4"),
      sorted_index_level_assignment = TRUE,
      exact_historical_observations_recoverable = FALSE
    )
  )
}

full_smoke_parameters <- function(parameters, generator, design) {
  parameters <- parameters
  if (identical(generator, "controlled_overlap")) {
    parameters$n <- as.integer(design$execution$smoke$n)
    parameters$p_numeric_noise <- min(as.integer(parameters$p_numeric_noise %||% 0L), 2L)
    parameters$q_categorical_noise <- min(as.integer(parameters$q_categorical_noise %||% 0L), 2L)
  } else if (identical(generator, "moons")) {
    parameters$n_per_moon <- as.integer(design$execution$smoke$n_per_moon)
    parameters$q_categorical_noise <- min(as.integer(parameters$q_categorical_noise %||% 0L), 1L)
  } else if (generator %in% c(
    "submitted_simulation_1",
    "submitted_simulation_2",
    "legacy_simulation_1",
    "legacy_simulation_2"
  )) {
    parameters$n <- as.integer(design$execution$smoke$n)
  }
  parameters
}

full_generate_scenario <- function(task, design, smoke = FALSE) {
  parameters <- task$parameters
  if (isTRUE(smoke)) {
    parameters <- full_smoke_parameters(parameters, task$generator, design)
  }
  generated <- switch(
    task$generator,
    controlled_overlap = full_make_controlled_overlap(parameters, task$data_seed),
    moons = full_make_moons(parameters, task$data_seed),
    submitted_simulation_1 = full_make_submitted_simulation_1(parameters, task$data_seed),
    submitted_simulation_2 = full_make_submitted_simulation_2(parameters, task$data_seed),
    legacy_simulation_1 = full_make_legacy_simulation_1(parameters, task$data_seed),
    legacy_simulation_2 = full_make_legacy_simulation_2(parameters, task$data_seed),
    stop("Unknown scenario generator: ", task$generator, call. = FALSE)
  )
  generated$parameters <- parameters
  generated
}
