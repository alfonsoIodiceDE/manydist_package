#!/usr/bin/env Rscript

script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
}

suite_dir <- dirname(script_path())
revision_dir <- dirname(suite_dir)
repo_root <- dirname(dirname(revision_dir))
full_dir <- file.path(revision_dir, "full")
design_path <- file.path(suite_dir, "design.yml")
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("The smoke summarizer requires `yaml`.", call. = FALSE)
}
design <- yaml::read_yaml(design_path)
design_tag <- substr(unname(tools::md5sum(design_path)), 1L, 10L)
input_dir <- file.path(suite_dir, "results", "smoke", design_tag)
current_source_files <- c(
  design = design_path,
  runner = file.path(suite_dir, "run_task.R"),
  generator = file.path(suite_dir, "R", "generator.R"),
  gated_distance = file.path(suite_dir, "R", "gated_distance.R"),
  bridge = file.path(revision_dir, "R", "manydist_bridge.R"),
  gate_helpers = file.path(revision_dir, "R", "distance_r1.R"),
  spectral_engine = file.path(full_dir, "R", "engine.R")
)
current_source_hashes <- unname(tools::md5sum(current_source_files))
names(current_source_hashes) <- names(current_source_files)

configurations <- names(design$generator$configurations)
states <- unlist(design$generator$interaction_states, use.names = FALSE)
replicates <- seq_len(as.integer(design$execution$smoke_replicates))
expected <- expand.grid(
  configuration = configurations,
  interaction_state = states,
  replicate = replicates,
  stringsAsFactors = FALSE
)
expected$file <- file.path(
  input_dir,
  expected$configuration,
  expected$interaction_state,
  sprintf("replicate-%02d.rds", expected$replicate)
)
missing <- expected$file[!file.exists(expected$file)]
if (length(missing)) {
  stop(
    "Missing ", length(missing), " expected artifact(s):\n",
    paste(missing, collapse = "\n"),
    call. = FALSE
  )
}

artifacts <- lapply(expected$file, readRDS)
valid <- vapply(seq_along(artifacts), function(index) {
  artifact <- artifacts[[index]]
  identical(artifact$analysis_id, as.character(design$analysis_id)) &&
    identical(artifact$design_tag, design_tag) &&
    identical(artifact$configuration, expected$configuration[[index]]) &&
    identical(artifact$interaction_state, expected$interaction_state[[index]]) &&
    identical(as.integer(artifact$replicate), expected$replicate[[index]]) &&
    identical(artifact$manydist_version, as.character(design$manydist$version)) &&
    identical(artifact$source_commit, as.character(design$manydist$source_commit)) &&
    identical(artifact$source_hashes, current_source_hashes)
}, logical(1))
if (!all(valid)) {
  stop("At least one artifact does not match the declared design.", call. = FALSE)
}

metrics <- do.call(rbind, lapply(artifacts, `[[`, "results"))
rownames(metrics) <- NULL
gates <- do.call(rbind, lapply(artifacts, function(artifact) {
  detected <- artifact$gate$detected
  data.frame(
    configuration = artifact$configuration,
    interaction_state = artifact$interaction_state,
    replicate = artifact$replicate,
    categorical_variable = names(detected),
    detected = unname(detected),
    p_value = unname(artifact$gate$p_value),
    statistic = unname(artifact$gate$statistic),
    active_count = artifact$gate$active_count,
    numeric_mean_distinct = artifact$distance_diagnostics$numeric_mean_distinct,
    added_mean_distinct = artifact$distance_diagnostics$added_mean_distinct,
    added_to_numeric_mean_ratio =
      artifact$distance_diagnostics$added_mean_distinct /
      artifact$distance_diagnostics$numeric_mean_distinct,
    stringsAsFactors = FALSE
  )
}))
rownames(gates) <- NULL

gated <- metrics[metrics$method == "ab_gated", c(
  "configuration", "interaction_state", "replicate", "affinity", "ARI"
)]
names(gated)[names(gated) == "ARI"] <- "ARI_gated"
baseline <- metrics[metrics$method == "ab_no_interaction", c(
  "configuration", "interaction_state", "replicate", "affinity", "ARI"
)]
names(baseline)[names(baseline) == "ARI"] <- "ARI_no_interaction"
contrasts <- merge(
  gated,
  baseline,
  by = c("configuration", "interaction_state", "replicate", "affinity"),
  sort = FALSE
)
contrasts$ARI_gain <- contrasts$ARI_gated - contrasts$ARI_no_interaction
active_counts <- unique(gates[c(
  "configuration", "interaction_state", "replicate", "active_count",
  "numeric_mean_distinct", "added_mean_distinct",
  "added_to_numeric_mean_ratio"
)])
contrasts <- merge(
  contrasts,
  active_counts,
  by = c("configuration", "interaction_state", "replicate"),
  sort = FALSE
)

configuration_order <- configurations
state_order <- states
affinity_order <- names(design$spectral$affinities)
method_order <- c(
  "ab_gated", "ab_no_interaction", "gower", "modified_gower",
  "standardized_euclidean_onehot"
)
metrics <- metrics[order(
  match(metrics$configuration, configuration_order),
  match(metrics$interaction_state, state_order),
  metrics$replicate,
  match(metrics$method, method_order),
  match(metrics$affinity, affinity_order)
), ]
gates <- gates[order(
  match(gates$configuration, configuration_order),
  match(gates$interaction_state, state_order),
  gates$replicate,
  gates$categorical_variable
), ]
contrasts <- contrasts[order(
  match(contrasts$configuration, configuration_order),
  match(contrasts$interaction_state, state_order),
  contrasts$replicate,
  match(contrasts$affinity, affinity_order)
), ]

write.csv(metrics, file.path(input_dir, "metrics.csv"), row.names = FALSE)
write.csv(gates, file.path(input_dir, "gates.csv"), row.names = FALSE)
write.csv(contrasts, file.path(input_dir, "gated_contrasts.csv"), row.names = FALSE)

cat("Validated and summarized ", length(artifacts), " artifact(s).\n", sep = "")
cat("Design tag: ", design_tag, "\n", sep = "")
cat("Output: ", input_dir, "\n", sep = "")
print(contrasts[c(
  "configuration", "interaction_state", "affinity", "ARI_gated",
  "ARI_no_interaction", "ARI_gain", "active_count"
)], row.names = FALSE)
