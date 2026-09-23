#!/usr/bin/env Rscript

# Post-hoc interaction-contribution diagnostics for the three-seed historical
# Simulation 1 pilot.  The diagnostic quantities use estimated partitions;
# truth is used only to retain the simulation ARI comparison.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "spectral_cluster_experiments/revision_r1/full/run_legacy_sim1_posthoc.R",
    mustWork = TRUE
  )
}

full_dir <- dirname(script_path)
revision_dir <- dirname(full_dir)
repo_root <- dirname(dirname(revision_dir))
output_dir <- file.path(
  full_dir,
  "results",
  "posthoc",
  "legacy_sim1_three_seed"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(full_dir, "R", "design.R"))
source(file.path(full_dir, "R", "scenarios.R"))

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

design <- full_read_design(file.path(full_dir, "design.yml"))
data_seeds <- c(26092201L, 26092301L, 26092401L)
parameters <- list(n = 500L, p_numeric = 15L, q_categorical = 15L)
gate_seed_offset <- as.integer(
  design$distance$categorical$interaction_gate$seed_offset
)
cluster_seed_offset <- as.integer(
  design$execution$cluster_seed_offsets$gaussian
)

alignment_statistics <- function(distance, partition) {
  distance <- full_validate_distance(distance)
  partition <- as.character(partition)
  upper <- upper.tri(distance)
  same_cluster <- outer(partition, partition, `==`)[upper]
  values <- distance[upper]
  within <- mean(values[same_cluster])
  between <- mean(values[!same_cluster])
  overall_sd <- stats::sd(values)
  data.frame(
    within_mean = within,
    between_mean = between,
    between_within_ratio = if (within > 0) between / within else NA_real_,
    between_within_contrast = between - within,
    standardized_contrast = if (overall_sd > 0) {
      (between - within) / overall_sd
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

cluster_distance <- function(distance, seed, k) {
  affinity <- full_affinity(
    distance,
    "gaussian",
    design$spectral$affinities$gaussian
  )
  fit <- full_spectral_cluster(
    affinity$matrix,
    k = k,
    seed = seed,
    kmeans_specification = design$spectral$kmeans
  )
  list(cluster = fit$cluster, sigma = affinity$scale)
}

mean_distance_matrices <- function(matrices, template) {
  if (!length(matrices)) {
    return(matrix(0, nrow(template), ncol(template), dimnames = dimnames(template)))
  }
  Reduce(`+`, matrices) / length(matrices)
}

replicate_rows <- list()
component_rows <- list()

for (replicate_index in seq_along(data_seeds)) {
  data_seed <- data_seeds[[replicate_index]]
  generated <- full_make_legacy_simulation_1(parameters, data_seed)
  truth <- droplevels(factor(generated$truth))

  gated_fit <- scmix_gated_active_distance(
    data = generated$x,
    gamma = as.numeric(design$distance$defaults$gamma),
    prop_nn = as.numeric(design$distance$defaults$prop_nn),
    score = as.character(design$distance$defaults$score),
    decision = as.character(design$distance$defaults$decision),
    gate_permutations = as.integer(
      design$distance$categorical$interaction_gate$permutations
    ),
    gate_alpha = as.numeric(
      design$distance$categorical$interaction_gate$alpha
    ),
    gate_seed = data_seed + gate_seed_offset
  )

  baseline_distance <- gated_fit$distance_no_interaction
  gated_distance <- gated_fit$distance
  active_interaction <- gated_fit$components$interaction_active_average
  active_names <- gated_fit$diagnostics$active_names
  cluster_seed <- data_seed + cluster_seed_offset
  k <- nlevels(truth)

  baseline_cluster <- cluster_distance(baseline_distance, cluster_seed, k)
  gated_cluster <- cluster_distance(gated_distance, cluster_seed, k)
  baseline_alignment <- alignment_statistics(
    active_interaction,
    baseline_cluster$cluster
  )
  gated_alignment <- alignment_statistics(
    active_interaction,
    gated_cluster$cluster
  )

  upper <- upper.tri(baseline_distance)
  redundancy_spearman <- suppressWarnings(stats::cor(
    baseline_distance[upper],
    active_interaction[upper],
    method = "spearman"
  ))

  replicate_components <- list()
  for (component_index in seq_along(active_names)) {
    variable <- active_names[[component_index]]
    held_out_interaction <- gated_fit$components$interaction_observation[[variable]]
    remaining_names <- setdiff(active_names, variable)
    remaining_average <- mean_distance_matrices(
      gated_fit$components$interaction_observation[remaining_names],
      baseline_distance
    )
    leave_one_out_distance <- full_validate_distance(
      baseline_distance +
        gated_fit$diagnostics$numeric_mean_distinct * remaining_average,
      paste("Leave-one-interaction-out distance for", variable)
    )
    leave_one_out_cluster <- cluster_distance(
      leave_one_out_distance,
      cluster_seed,
      k
    )
    held_out_alignment <- alignment_statistics(
      held_out_interaction,
      leave_one_out_cluster$cluster
    )
    direct_baseline_alignment <- alignment_statistics(
      held_out_interaction,
      baseline_cluster$cluster
    )
    direct_gated_alignment <- alignment_statistics(
      held_out_interaction,
      gated_cluster$cluster
    )

    replicate_components[[component_index]] <- data.frame(
      replicate = replicate_index,
      data_seed = data_seed,
      variable = variable,
      gate_p_value = gated_fit$gate$p_value[[variable]],
      gate_statistic = gated_fit$gate$statistic[[variable]],
      baseline_alignment_ratio = direct_baseline_alignment$between_within_ratio,
      gated_alignment_ratio = direct_gated_alignment$between_within_ratio,
      leave_one_out_alignment_ratio = held_out_alignment$between_within_ratio,
      leave_one_out_standardized_contrast = held_out_alignment$standardized_contrast,
      leave_one_out_partition_ARI_to_gated = mclust::adjustedRandIndex(
        leave_one_out_cluster$cluster,
        gated_cluster$cluster
      ),
      leave_one_out_ARI_to_truth = mclust::adjustedRandIndex(
        leave_one_out_cluster$cluster,
        truth
      ),
      stringsAsFactors = FALSE
    )
  }
  replicate_components <- do.call(rbind, replicate_components)
  component_rows[[replicate_index]] <- replicate_components

  replicate_rows[[replicate_index]] <- data.frame(
    replicate = replicate_index,
    data_seed = data_seed,
    active_count = length(active_names),
    baseline_ARI_to_truth = mclust::adjustedRandIndex(
      baseline_cluster$cluster,
      truth
    ),
    gated_ARI_to_truth = mclust::adjustedRandIndex(
      gated_cluster$cluster,
      truth
    ),
    gated_minus_baseline_ARI = mclust::adjustedRandIndex(
      gated_cluster$cluster,
      truth
    ) - mclust::adjustedRandIndex(baseline_cluster$cluster, truth),
    baseline_gated_partition_ARI = mclust::adjustedRandIndex(
      baseline_cluster$cluster,
      gated_cluster$cluster
    ),
    active_interaction_alignment_ratio_baseline =
      baseline_alignment$between_within_ratio,
    active_interaction_alignment_ratio_gated =
      gated_alignment$between_within_ratio,
    active_interaction_standardized_contrast_baseline =
      baseline_alignment$standardized_contrast,
    active_interaction_standardized_contrast_gated =
      gated_alignment$standardized_contrast,
    interaction_baseline_spearman = redundancy_spearman,
    mean_leave_one_out_alignment_ratio = mean(
      replicate_components$leave_one_out_alignment_ratio
    ),
    min_leave_one_out_alignment_ratio = min(
      replicate_components$leave_one_out_alignment_ratio
    ),
    proportion_leave_one_out_aligned = mean(
      replicate_components$leave_one_out_alignment_ratio > 1
    ),
    mean_leave_one_out_partition_ARI_to_gated = mean(
      replicate_components$leave_one_out_partition_ARI_to_gated
    ),
    baseline_sigma = baseline_cluster$sigma,
    gated_sigma = gated_cluster$sigma,
    stringsAsFactors = FALSE
  )
}

replicate_results <- do.call(rbind, replicate_rows)
component_results <- do.call(rbind, component_rows)

numeric_columns <- names(replicate_results)[vapply(
  replicate_results,
  is.numeric,
  logical(1)
)]
numeric_columns <- setdiff(numeric_columns, c("replicate", "data_seed"))
summary_results <- do.call(rbind, lapply(numeric_columns, function(variable) {
  values <- replicate_results[[variable]]
  data.frame(
    quantity = variable,
    mean = mean(values),
    sd = stats::sd(values),
    min = min(values),
    max = max(values),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(
  replicate_results,
  file.path(output_dir, "replicate_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  component_results,
  file.path(output_dir, "component_leave_one_out.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_results,
  file.path(output_dir, "summary.csv"),
  row.names = FALSE
)

cat("Wrote post-hoc diagnostics to:\n", output_dir, "\n", sep = "")
print(replicate_results, digits = 4)
