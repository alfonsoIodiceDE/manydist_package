structured_blocks_correlation_from_loadings <- function(loadings) {
  loadings <- as.matrix(loadings)
  communalities <- rowSums(loadings^2)
  if (any(communalities >= 1)) {
    stop("Every loading row must have squared norm below one.", call. = FALSE)
  }
  covariance <- tcrossprod(loadings) + diag(1 - communalities)
  covariance <- stats::cov2cor(covariance)
  if (min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Constructed correlation matrix is not positive definite.", call. = FALSE)
  }
  covariance
}

structured_blocks_loading_matrices <- function(p, loading) {
  if (p %% 3L != 0L) {
    stop("The pilot requires `p` divisible by three.", call. = FALSE)
  }
  block_size <- p / 3L

  sequential <- matrix(0, p, 3L)
  sequential[cbind(seq_len(p), rep(seq_len(3L), each = block_size))] <- loading

  interleaved <- matrix(0, p, 3L)
  interleaved[cbind(seq_len(p), rep(seq_len(3L), length.out = p))] <- loading

  signed <- matrix(0, p, 2L)
  signed[seq_len(p / 2L), 1L] <- loading
  signed[seq.int(p / 2L + 1L, p), 2L] <- loading
  signed[seq(2L, p / 2L, by = 2L), 1L] <- -loading
  signed[seq.int(p / 2L + 2L, p, by = 2L), 2L] <- -loading

  list(sequential = sequential, interleaved = interleaved, signed = signed)
}

structured_blocks_numeric_means <- function(p, scale) {
  list(
    rep(0, p),
    scale * rep(c(1, -0.25, 0.25), length.out = p),
    scale * rep(c(-0.5, 1, -1), length.out = p)
  )
}

structured_blocks_make_categorical <- function(
    regime,
    levels_per_variable,
    correlation_matrices,
    seed
) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("The structured-block generator requires `mvtnorm`.", call. = FALSE)
  }
  regime <- as.integer(regime)
  q <- length(levels_per_variable)
  latent <- matrix(NA_real_, length(regime), q)

  set.seed(seed)
  for (group in seq_along(correlation_matrices)) {
    rows <- which(regime == group)
    latent[rows, ] <- mvtnorm::rmvnorm(
      length(rows),
      mean = rep(0, q),
      sigma = correlation_matrices[[group]],
      method = "chol"
    )
  }

  # A fixed arbitrary relabelling removes any operational use of category
  # order.  The same relabelling is used in every association regime.
  set.seed(seed + 1L)
  label_permutations <- lapply(levels_per_variable, function(level_count) {
    sample.int(level_count)
  })
  categorical <- lapply(seq_len(q), function(variable) {
    level_count <- levels_per_variable[[variable]]
    boundaries <- c(
      -Inf,
      stats::qnorm(seq_len(level_count - 1L) / level_count),
      Inf
    )
    original_code <- cut(
      latent[, variable],
      breaks = boundaries,
      labels = FALSE,
      include.lowest = TRUE
    )
    permuted_code <- label_permutations[[variable]][original_code]
    factor(
      paste0("L", permuted_code),
      levels = paste0("L", seq_len(level_count))
    )
  })
  names(categorical) <- paste0("C", seq_len(q))
  as.data.frame(categorical, check.names = FALSE)
}

structured_blocks_nominal_association_patterns <- function(variable_indices) {
  variable_count <- length(variable_indices)
  if (!variable_count %in% c(4L, 8L)) {
    stop(
      paste(
        "The nominal-block generator requires exactly four or eight",
        "variables of each arity."
      ),
      call. = FALSE
    )
  }
  i <- variable_indices

  # With four variables, the three regimes are the three possible balanced
  # two-block partitions. This is the smallest design that retains three
  # genuinely distinct nominal association structures.
  if (variable_count == 4L) {
    return(list(
      list(i[c(1, 2)], i[c(3, 4)]),
      list(i[c(1, 3)], i[c(2, 4)]),
      list(i[c(1, 4)], i[c(2, 3)])
    ))
  }

  # Preserve the original eight-variable construction for reproducibility of
  # the pilot artifacts.
  list(
    list(i[c(1, 2, 3, 4)], i[c(5, 6, 7, 8)]),
    list(i[c(1, 3, 5, 7)], i[c(2, 4, 6, 8)]),
    list(i[c(1, 2, 5, 6)], i[c(3, 4, 7, 8)])
  )
}

structured_blocks_make_nominal_categorical <- function(
    regime,
    levels_per_variable,
    copy_probability,
    profile_strength,
    seed
) {
  if (length(copy_probability) != 1L || !is.finite(copy_probability) ||
      copy_probability < 0 || copy_probability > 1) {
    stop("`copy_probability` must lie in [0, 1].", call. = FALSE)
  }
  if (length(profile_strength) != 1L || !is.finite(profile_strength) ||
      profile_strength < 0 || profile_strength >= 1) {
    stop("`profile_strength` must lie in [0, 1).", call. = FALSE)
  }
  regime <- as.integer(regime)
  q <- length(levels_per_variable)
  if (!identical(sort(unique(levels_per_variable)), 2:4)) {
    stop("The pilot expects categorical arities two, three, and four.",
         call. = FALSE)
  }

  set.seed(seed)
  label_permutations <- lapply(levels_per_variable, sample.int)
  values <- matrix(NA_integer_, nrow = length(regime), ncol = q)
  patterns_by_arity <- lapply(2:4, function(level_count) {
    structured_blocks_nominal_association_patterns(
      which(levels_per_variable == level_count)
    )
  })
  names(patterns_by_arity) <- as.character(2:4)

  for (row in seq_along(regime)) {
    association_regime <- regime[[row]]
    for (level_count in 2:4) {
      blocks <- patterns_by_arity[[as.character(level_count)]][[association_regime]]
      profile <- rep((1 - profile_strength) / level_count, level_count)
      favored_level <- (association_regime - 1L) %% level_count + 1L
      profile[[favored_level]] <- profile[[favored_level]] + profile_strength
      for (block in blocks) {
        shared_code <- sample.int(level_count, 1L, prob = profile)
        copied <- stats::runif(length(block)) < copy_probability
        raw_codes <- ifelse(
          copied,
          shared_code,
          sample.int(
            level_count,
            length(block),
            replace = TRUE,
            prob = profile
          )
        )
        for (position in seq_along(block)) {
          variable <- block[[position]]
          values[row, variable] <- label_permutations[[variable]][raw_codes[[position]]]
        }
      }
    }
  }

  categorical <- lapply(seq_len(q), function(variable) {
    level_count <- levels_per_variable[[variable]]
    factor(
      paste0("L", values[, variable]),
      levels = paste0("L", seq_len(level_count))
    )
  })
  names(categorical) <- paste0("C", seq_len(q))
  as.data.frame(categorical, check.names = FALSE)
}

structured_blocks_make_pair <- function(specification, seed) {
  n_per_group <- as.integer(specification$n_per_group)
  p <- as.integer(specification$p_numeric)
  categorical_levels <- as.integer(unlist(specification$categorical_levels))
  q <- length(categorical_levels)
  if (n_per_group %% 3L != 0L) {
    stop("`n_per_group` must be divisible by three for the exact decoupled arm.",
         call. = FALSE)
  }
  if (p %% 3L != 0L || q %% 3L != 0L) {
    stop("Numerical and categorical dimensions must be divisible by three.",
         call. = FALSE)
  }

  group <- rep(seq_len(3L), each = n_per_group)
  n <- length(group)

  numeric_loadings <- structured_blocks_loading_matrices(
    p,
    as.numeric(specification$numeric_loading)
  )
  numeric_correlations <- lapply(
    numeric_loadings,
    structured_blocks_correlation_from_loadings
  )
  numeric_means <- structured_blocks_numeric_means(
    p,
    as.numeric(specification$numeric_mean_scale)
  )

  set.seed(seed)
  numeric <- matrix(NA_real_, n, p)
  for (g in seq_len(3L)) {
    rows <- which(group == g)
    numeric[rows, ] <- mvtnorm::rmvnorm(
      length(rows),
      mean = numeric_means[[g]],
      sigma = numeric_correlations[[g]],
      method = "chol"
    )
  }
  colnames(numeric) <- paste0("V", seq_len(p))
  numeric <- as.data.frame(numeric, check.names = FALSE)

  categorical_regime <- rep(seq_len(3L), each = n_per_group)
  categorical_model <- as.character(specification$categorical_model)
  if (identical(categorical_model, "nominal_association_blocks")) {
    categorical_bank <- structured_blocks_make_nominal_categorical(
      categorical_regime,
      categorical_levels,
      copy_probability = as.numeric(specification$categorical_copy_probability),
      profile_strength = as.numeric(specification$categorical_profile_strength),
      seed = seed + 1000L
    )
  } else if (identical(categorical_model, "gaussian_copula")) {
    categorical_loadings <- structured_blocks_loading_matrices(
      q,
      as.numeric(specification$categorical_loading)
    )
    categorical_correlations <- lapply(
      categorical_loadings,
      structured_blocks_correlation_from_loadings
    )
    categorical_bank <- structured_blocks_make_categorical(
      categorical_regime,
      categorical_levels,
      categorical_correlations,
      seed = seed + 1000L
    )
  } else {
    stop("Unknown categorical model: ", categorical_model, call. = FALSE)
  }

  # The aligned arm matches association regime h to numerical group g = h.
  aligned_index <- seq_len(n)

  # The decoupled arm uses the same categorical rows exactly once, while each
  # numerical group receives an equal number from every association regime.
  set.seed(seed + 2000L)
  regime_rows <- lapply(seq_len(3L), function(h) {
    sample(which(categorical_regime == h))
  })
  chunk_size <- n_per_group / 3L
  decoupled_index <- integer(n)
  for (g in seq_len(3L)) {
    target_rows <- which(group == g)
    selected <- unlist(lapply(seq_len(3L), function(h) {
      start <- (g - 1L) * chunk_size + 1L
      stop <- g * chunk_size
      regime_rows[[h]][start:stop]
    }), use.names = FALSE)
    decoupled_index[target_rows] <- sample(selected)
  }

  make_arm <- function(index, arm) {
    regime <- categorical_regime[index]
    categorical_arm <- categorical_bank[index, , drop = FALSE]
    rownames(categorical_arm) <- NULL
    combined <- data.frame(numeric, categorical_arm, check.names = FALSE)
    rownames(combined) <- NULL
    list(
      x = combined,
      truth = factor(group, levels = seq_len(3L), labels = paste0("G", seq_len(3L))),
      categorical_regime = factor(
        regime,
        levels = seq_len(3L),
        labels = paste0("H", seq_len(3L))
      ),
      arm = arm
    )
  }

  list(
    aligned = make_arm(aligned_index, "aligned"),
    decoupled = make_arm(decoupled_index, "decoupled"),
    diagnostics = list(
      numeric_correlations = numeric_correlations,
      categorical_model = categorical_model,
      categorical_levels = categorical_levels,
      numeric_means = numeric_means
    )
  )
}
