# Paper-specific gated extension of a pure manydist custom baseline.
#
# Source revision_r1/R/manydist_bridge.R and revision_r1/R/distance_r1.R
# before this file. The package computes the complete no-interaction baseline
# and numerical geometry. Its internal delta_int_knn() computes category-level
# interaction deltas. This file only applies the paper-specific multiplicity
# gate and active-component aggregation.

moon_scaling_gated_distance <- function(
    data,
    gamma = 1,
    prop_nn = 0.10,
    score = "ba",
    decision = "prior_corrected",
    gate_permutations = 99L,
    gate_alpha = 0.01,
    gate_seed = 1L
) {
  required_helpers <- c(
    ".scmix_nn_order",
    ".scmix_interaction_delta",
    ".scmix_permutation_max_gate",
    ".scmix_observation_distance",
    ".scmix_mean_distinct"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function"
  )]
  if (length(missing_helpers)) {
    stop(
      "Source the revision bridge and distance helpers first; missing: ",
      paste(missing_helpers, collapse = ", "),
      call. = FALSE
    )
  }
  if (!requireNamespace("manydist", quietly = TRUE)) {
    stop("The `manydist` package is required.", call. = FALSE)
  }

  data <- as.data.frame(data)
  numeric_names <- names(data)[vapply(data, is.numeric, logical(1))]
  categorical_names <- names(data)[vapply(
    data,
    function(x) is.factor(x) || is.character(x),
    logical(1)
  )]
  if (!length(numeric_names) || !length(categorical_names) ||
      length(c(numeric_names, categorical_names)) != ncol(data)) {
    stop("Data must contain only numerical and nominal categorical predictors.",
         call. = FALSE)
  }
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0 || gamma > 1) {
    stop("`gamma` must be a finite scalar in [0, 1].", call. = FALSE)
  }
  if (length(prop_nn) != 1L || !is.finite(prop_nn) || prop_nn <= 0 || prop_nn > 1) {
    stop("`prop_nn` must lie in (0, 1].", call. = FALSE)
  }
  score <- match.arg(score, c("ba", "logloss"))
  decision <- match.arg(decision, c("prior_corrected", "posterior"))
  if (length(gate_alpha) != 1L || !is.finite(gate_alpha) ||
      gate_alpha <= 0 || gate_alpha >= 1) {
    stop("`gate_alpha` must lie in (0, 1).", call. = FALSE)
  }
  if (length(gate_permutations) != 1L || !is.finite(gate_permutations) ||
      gate_permutations != as.integer(gate_permutations) || gate_permutations < 1L) {
    stop("`gate_permutations` must be a positive integer.", call. = FALSE)
  }
  gate_permutations <- as.integer(gate_permutations)
  if (1 / (gate_permutations + 1) > gate_alpha) {
    stop("The requested permutation count cannot attain `gate_alpha`.", call. = FALSE)
  }
  if (length(gate_seed) != 1L || !is.finite(gate_seed) ||
      gate_seed != as.integer(gate_seed) || gate_seed < 0) {
    stop("`gate_seed` must be a nonnegative integer.", call. = FALSE)
  }
  gate_seed <- as.integer(gate_seed)

  ncomp <- length(numeric_names)
  package_arguments <- list(
    preset = "custom",
    method_num = "pc_scores",
    method_cat = "tvd",
    commensurable = FALSE,
    interaction = FALSE,
    ncomp = ncomp
  )

  baseline_started <- proc.time()[[3L]]
  baseline_fit <- do.call(
    manydist::mdist,
    c(list(x = data), package_arguments)
  )
  baseline <- .scmix_validate_distance(
    as.matrix(baseline_fit$distance),
    "Pure manydist no-interaction baseline"
  )
  baseline_seconds <- proc.time()[[3L]] - baseline_started

  numeric_started <- proc.time()[[3L]]
  numeric_fit <- do.call(
    manydist::mdist,
    c(list(x = data[numeric_names]), package_arguments)
  )
  numeric_distance <- .scmix_validate_distance(
    as.matrix(numeric_fit$distance),
    "Pure manydist numerical distance"
  )
  numeric_seconds <- proc.time()[[3L]] - numeric_started

  categorical_started <- proc.time()[[3L]]
  categorical_fit <- do.call(
    manydist::mdist,
    c(list(x = data[categorical_names]), package_arguments)
  )
  categorical_distance <- .scmix_validate_distance(
    as.matrix(categorical_fit$distance),
    "Pure manydist categorical distance"
  )
  categorical_seconds <- proc.time()[[3L]] - categorical_started
  decomposition_error <- max(abs(
    baseline - (numeric_distance + categorical_distance)
  ))
  if (decomposition_error > 1e-8) {
    stop("The package baseline does not equal its numerical and categorical calls.",
         call. = FALSE)
  }

  categorical <- lapply(data[categorical_names], function(x) droplevels(factor(x)))
  categorical <- as.data.frame(categorical, stringsAsFactors = TRUE)
  nn_order <- .scmix_nn_order(numeric_distance)
  interaction_level_delta <- lapply(categorical, function(labels) {
    .scmix_interaction_delta(
      numeric_distance = numeric_distance,
      labels = labels,
      prop_nn = prop_nn,
      score = score,
      decision = decision,
      nn_order = nn_order
    )
  })
  names(interaction_level_delta) <- categorical_names

  gate_started <- proc.time()[[3L]]
  gate <- .scmix_permutation_max_gate(
    numeric_distance = numeric_distance,
    categorical = categorical,
    observed_level_deltas = interaction_level_delta,
    prop_nn = prop_nn,
    score = score,
    decision = decision,
    nn_order = nn_order,
    permutations = gate_permutations,
    alpha = gate_alpha,
    seed = gate_seed
  )
  gate_seconds <- proc.time()[[3L]] - gate_started

  interaction_observation <- lapply(seq_along(categorical), function(j) {
    .scmix_observation_distance(interaction_level_delta[[j]], categorical[[j]])
  })
  names(interaction_observation) <- categorical_names
  active_names <- names(gate$detected)[gate$detected]
  active_count <- length(active_names)
  zero <- matrix(0, nrow(data), nrow(data), dimnames = dimnames(baseline))
  active_average <- if (!active_count) {
    zero
  } else {
    Reduce(`+`, interaction_observation[active_names], init = zero) / active_count
  }
  active_average <- .scmix_validate_distance(
    active_average,
    "Active-average interaction distance"
  )

  numeric_mean <- .scmix_mean_distinct(numeric_distance)
  interaction_weight <- gamma * numeric_mean
  added_interaction <- interaction_weight * active_average
  gated <- if (!active_count || gamma == 0) {
    baseline
  } else {
    .scmix_validate_distance(
      baseline + added_interaction,
      "Gated interaction distance"
    )
  }
  if (!active_count && !identical(gated, baseline)) {
    stop("No active gate must recover the package baseline exactly.", call. = FALSE)
  }

  structure(
    list(
      distance = gated,
      distance_no_interaction = baseline,
      components = list(
        numeric = numeric_distance,
        categorical = categorical_distance,
        interaction_level_delta = interaction_level_delta,
        interaction_observation = interaction_observation,
        interaction_active_average = active_average,
        interaction_added = added_interaction
      ),
      gate = gate,
      diagnostics = list(
        n = nrow(data),
        p_numeric = length(numeric_names),
        q_categorical = length(categorical_names),
        active_names = active_names,
        active_count = active_count,
        gamma = gamma,
        numeric_mean_distinct = numeric_mean,
        categorical_mean_distinct = .scmix_mean_distinct(categorical_distance),
        active_average_mean_distinct = .scmix_mean_distinct(active_average),
        added_mean_distinct = .scmix_mean_distinct(added_interaction),
        final_mean_distinct = .scmix_mean_distinct(gated),
        decomposition_error = decomposition_error,
        baseline_seconds = baseline_seconds,
        numeric_seconds = numeric_seconds,
        categorical_seconds = categorical_seconds,
        gate_seconds = gate_seconds
      ),
      configuration = list(
        baseline = package_arguments,
        gamma = gamma,
        aggregation = "mean over multiplicity-gated active interaction components",
        prop_nn = prop_nn,
        score = score,
        decision = decision,
        gate_permutations = gate_permutations,
        gate_alpha = gate_alpha,
        gate_seed = gate_seed
      )
    ),
    class = "moon_scaling_gated_distance"
  )
}

# The implementation is data-agnostic despite its historical location in the
# Simulation 3 folder.  The full revision runner uses this explicit general
# name so every scenario shares exactly the construction validated here.
scmix_gated_active_distance <- moon_scaling_gated_distance
