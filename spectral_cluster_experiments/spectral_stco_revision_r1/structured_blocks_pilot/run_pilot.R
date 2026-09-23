#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "spectral_cluster_experiments/revision_r1/structured_blocks_pilot/run_pilot.R",
    mustWork = TRUE
  )
}

pilot_dir <- dirname(script_path)
revision_dir <- dirname(pilot_dir)
full_dir <- file.path(revision_dir, "full")
repo_root <- dirname(dirname(revision_dir))
output_dir <- file.path(pilot_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("devtools", "mclust", "mvtnorm", "yaml")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)
source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))
source(file.path(revision_dir, "moon_scaling", "R", "gated_distance.R"))
source(file.path(full_dir, "R", "engine.R"))
source(file.path(pilot_dir, "R", "generator.R"))

specification <- yaml::read_yaml(file.path(pilot_dir, "design.yml"))
data_specification <- specification$data
distance_specification <- specification$distance
spectral_specification <- specification$spectral
data_seeds <- as.integer(unlist(data_specification$seeds))
arms <- as.character(unlist(data_specification$arms))

cramers_v <- function(x, y) {
  observed <- table(x, y)
  n <- sum(observed)
  dimension <- min(nrow(observed) - 1L, ncol(observed) - 1L)
  if (!n || dimension <= 0) return(0)
  expected <- outer(rowSums(observed), colSums(observed)) / n
  valid <- expected > 0
  statistic <- sum((observed[valid] - expected[valid])^2 / expected[valid])
  sqrt((statistic / n) / dimension)
}

categorical_association_matrix <- function(data) {
  q <- ncol(data)
  result <- diag(1, q)
  for (i in seq_len(q - 1L)) {
    for (j in seq.int(i + 1L, q)) {
      result[i, j] <- result[j, i] <- cramers_v(data[[i]], data[[j]])
    }
  }
  result
}

mean_pairwise_matrix_distance <- function(matrices) {
  pairs <- utils::combn(seq_along(matrices), 2L)
  mean(apply(pairs, 2L, function(pair) {
    sqrt(sum((matrices[[pair[[1L]]]] - matrices[[pair[[2L]]]])^2))
  }))
}

categorical_structure_diagnostics <- function(generated) {
  categorical <- generated$x[vapply(generated$x, is.factor, logical(1))]
  truth <- generated$truth
  association <- lapply(levels(truth), function(level) {
    categorical_association_matrix(categorical[truth == level, , drop = FALSE])
  })

  marginal_deviations <- unlist(lapply(categorical, function(variable) {
    overall <- prop.table(table(variable))
    unlist(lapply(levels(truth), function(level) {
      within <- prop.table(table(variable[truth == level]))
      max(abs(within[names(overall)] - overall))
    }))
  }))
  data.frame(
    association_matrix_separation = mean_pairwise_matrix_distance(association),
    mean_marginal_deviation = mean(marginal_deviations),
    max_marginal_deviation = max(marginal_deviations),
    stringsAsFactors = FALSE
  )
}

cluster_distance <- function(distance, truth, seed) {
  affinity <- full_affinity(
    distance,
    "gaussian",
    list()
  )
  fit <- full_spectral_cluster(
    affinity$matrix,
    k = nlevels(truth),
    seed = seed,
    kmeans_specification = spectral_specification$kmeans
  )
  list(
    cluster = fit$cluster,
    ARI = mclust::adjustedRandIndex(fit$cluster, truth),
    sigma = affinity$scale
  )
}

alignment_ratio <- function(distance, partition) {
  distance <- full_validate_distance(distance)
  upper <- upper.tri(distance)
  same <- outer(partition, partition, `==`)[upper]
  values <- distance[upper]
  within <- mean(values[same])
  between <- mean(values[!same])
  if (within > 0) between / within else NA_real_
}

method_rows <- list()
gate_rows <- list()
structure_rows <- list()
row_index <- 0L
gate_index <- 0L
structure_index <- 0L

for (replicate_index in seq_along(data_seeds)) {
  data_seed <- data_seeds[[replicate_index]]
  generated_pair <- structured_blocks_make_pair(data_specification, data_seed)

  for (arm in arms) {
    generated <- generated_pair[[arm]]
    truth <- generated$truth
    data <- generated$x
    categorical <- data[vapply(data, is.factor, logical(1))]
    numeric <- data[vapply(data, is.numeric, logical(1))]
    cluster_seed <- data_seed + as.integer(spectral_specification$cluster_seed_offset)

    structure_index <- structure_index + 1L
    structure_rows[[structure_index]] <- cbind(
      data.frame(
        replicate = replicate_index,
        data_seed = data_seed,
        arm = arm,
        regime_truth_ARI = mclust::adjustedRandIndex(
          generated$categorical_regime,
          truth
        ),
        stringsAsFactors = FALSE
      ),
      categorical_structure_diagnostics(generated)
    )

    gated <- scmix_gated_active_distance(
      data = data,
      gamma = as.numeric(distance_specification$gamma),
      prop_nn = as.numeric(distance_specification$prop_nn),
      score = as.character(distance_specification$score),
      decision = as.character(distance_specification$decision),
      gate_permutations = as.integer(distance_specification$gate$permutations),
      gate_alpha = as.numeric(distance_specification$gate$alpha),
      gate_seed = data_seed + as.integer(distance_specification$gate$seed_offset)
    )

    distances <- list(
      Gated_active_average_SC = gated$distance,
      Package_no_interaction_SC = gated$distance_no_interaction,
      Gower_SC = as.matrix(manydist::mdist(data, preset = "gower")$distance),
      Modified_Gower_SC = as.matrix(manydist::mdist(data, preset = "mod_gower")$distance),
      Euclidean_onehot_SC = as.matrix(manydist::mdist(data, preset = "euclidean")$distance),
      Numeric_only_PC_SC = gated$components$numeric,
      Categorical_only_TVD_SC = gated$components$categorical
    )

    fitted_partitions <- list()
    for (method in names(distances)) {
      fitted <- cluster_distance(distances[[method]], truth, cluster_seed)
      fitted_partitions[[method]] <- fitted$cluster
      row_index <- row_index + 1L
      method_rows[[row_index]] <- data.frame(
        replicate = replicate_index,
        data_seed = data_seed,
        arm = arm,
        method = method,
        ARI = fitted$ARI,
        sigma = fitted$sigma,
        stringsAsFactors = FALSE
      )
    }

    if (!gated$diagnostics$active_count &&
        !identical(gated$distance, gated$distance_no_interaction)) {
      stop("A closed gate failed to reproduce the package baseline exactly.")
    }

    gate_index <- gate_index + 1L
    gate_rows[[gate_index]] <- data.frame(
      replicate = replicate_index,
      data_seed = data_seed,
      arm = arm,
      active_count = gated$diagnostics$active_count,
      active_names = paste(gated$diagnostics$active_names, collapse = ";"),
      minimum_adjusted_p = min(gated$gate$p_value),
      exact_baseline_recovery = identical(
        gated$distance,
        gated$distance_no_interaction
      ),
      baseline_gated_partition_ARI = mclust::adjustedRandIndex(
        fitted_partitions$Package_no_interaction_SC,
        fitted_partitions$Gated_active_average_SC
      ),
      interaction_alignment_ratio = if (gated$diagnostics$active_count) {
        alignment_ratio(
          gated$components$interaction_active_average,
          fitted_partitions$Gated_active_average_SC
        )
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
}

method_results <- do.call(rbind, method_rows)
gate_results <- do.call(rbind, gate_rows)
structure_results <- do.call(rbind, structure_rows)

method_summary <- stats::aggregate(
  ARI ~ arm + method,
  method_results,
  function(values) c(
    mean = mean(values),
    sd = stats::sd(values),
    min = min(values),
    max = max(values)
  )
)
method_summary <- data.frame(
  arm = method_summary$arm,
  method = method_summary$method,
  mean_ARI = method_summary$ARI[, "mean"],
  sd_ARI = method_summary$ARI[, "sd"],
  min_ARI = method_summary$ARI[, "min"],
  max_ARI = method_summary$ARI[, "max"],
  row.names = NULL,
  check.names = FALSE
)

utils::write.csv(method_results, file.path(output_dir, "method_results.csv"), row.names = FALSE)
utils::write.csv(method_summary, file.path(output_dir, "method_summary.csv"), row.names = FALSE)
utils::write.csv(gate_results, file.path(output_dir, "gate_results.csv"), row.names = FALSE)
utils::write.csv(structure_results, file.path(output_dir, "structure_results.csv"), row.names = FALSE)

cat("Method summary\n")
print(method_summary, digits = 4, row.names = FALSE)
cat("\nGate results\n")
print(gate_results, digits = 4, row.names = FALSE)
cat("\nStructure checks\n")
print(structure_results, digits = 4, row.names = FALSE)

