#!/usr/bin/env Rscript

# Exploratory, paired Simulation 2 check only. This does not alter the full
# revision manifest or establish the identity of the archived 50-run driver.

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE)
pilot_dir <- dirname(script)
revision_dir <- dirname(pilot_dir)
repo_root <- dirname(dirname(revision_dir))
full_dir <- file.path(revision_dir, "full")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript score_scale_aware_gate_smoke.R PATH_TO_PILOT_RDS", call. = FALSE)
}

artifact_path <- normalizePath(arguments[[1L]], mustWork = TRUE)
results_root <- paste0(normalizePath(file.path(pilot_dir, "results"), mustWork = TRUE),
                       .Platform$file.sep)
if (!startsWith(artifact_path, results_root) || !grepl("[.]rds$", artifact_path)) {
  stop("Input must be a saved paired pilot RDS under results/.", call. = FALSE)
}
relative_path <- substring(artifact_path, nchar(results_root) + 1L)
output_path <- file.path(pilot_dir, "scale_aware_gate_smoke", relative_path)
if (file.exists(output_path)) {
  cat("Exploratory sidecar already exists; leaving it unchanged: ", output_path, "\n", sep = "")
  quit(status = 0L)
}

source(file.path(full_dir, "R", "design.R"))
source(file.path(full_dir, "R", "engine.R"))
source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))
for (package in c("devtools", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("This scorer requires ", package, ".", call. = FALSE)
  }
}
devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)

design <- full_read_design(file.path(full_dir, "design.yml"))
artifact <- readRDS(artifact_path)
if (!identical(artifact$analysis_id, "sim2-code-reconstruction-paired-pilot") ||
    !identical(artifact$design_config_md5, attr(design, "config_md5")) ||
    nrow(artifact$data) != 500L) {
  stop("Input is not a paired n=500 Simulation 2 reconstruction artifact.", call. = FALSE)
}

data <- artifact$data
numeric_names <- names(data)[vapply(data, is.numeric, logical(1))]
categorical_names <- names(data)[vapply(data, is.factor, logical(1))]
stopifnot(length(numeric_names) == 6L, length(categorical_names) == 4L)
settings <- list(
  preset = "custom", method_num = "pc_scores", method_cat = "tvd",
  commensurable = FALSE, interaction = FALSE,
  ncomp = length(numeric_names)
)

# The no-interaction arm is the package result, not a reimplementation.
base <- full_validate_distance(as.matrix(do.call(manydist::mdist,
  c(list(x = data), settings))$distance), "Package custom baseline")
numeric_distance <- full_validate_distance(as.matrix(do.call(manydist::mdist,
  c(list(x = data[numeric_names]), settings))$distance), "Package numerical block")
categorical_distance <- full_validate_distance(as.matrix(do.call(manydist::mdist,
  c(list(x = data[categorical_names]), settings))$distance), "Package TVD block")
decomposition_error <- max(abs(base - numeric_distance - categorical_distance))
if (decomposition_error > 1e-8) {
  stop("The package baseline does not equal numerical plus TVD blocks.", call. = FALSE)
}

prop_nn <- as.numeric(design$distance$defaults$prop_nn)
score <- as.character(design$distance$defaults$score)
decision <- as.character(design$distance$defaults$decision)
gate_spec <- design$distance$categorical$interaction_gate
categorical <- lapply(data[categorical_names], function(x) droplevels(factor(x)))
nn_order <- .scmix_nn_order(numeric_distance)
delta <- lapply(categorical, function(labels) .scmix_interaction_delta(
  numeric_distance = numeric_distance, labels = labels, prop_nn = prop_nn,
  score = score, decision = decision, nn_order = nn_order
))

gate_started <- proc.time()[[3L]]
gate <- .scmix_permutation_max_gate(
  numeric_distance = numeric_distance,
  categorical = categorical,
  observed_level_deltas = delta,
  prop_nn = prop_nn,
  score = score,
  decision = decision,
  nn_order = nn_order,
  permutations = as.integer(gate_spec$permutations),
  alpha = as.numeric(gate_spec$alpha),
  seed = artifact$data_seed + as.integer(gate_spec$seed_offset)
)
gate_seconds <- proc.time()[[3L]] - gate_started

# This is Z_j %*% Delta_int,j %*% t(Z_j) at the category level, gated before
# projection. Gamma is exploratory; neither value is a selected primary method.
gated_interaction <- matrix(0, nrow(data), nrow(data))
for (name in categorical_names) {
  if (!isTRUE(gate$detected[[name]])) next
  z <- diag(nlevels(categorical[[name]]))[as.integer(categorical[[name]]), , drop = FALSE]
  gated_interaction <- gated_interaction + z %*% delta[[name]] %*% t(z)
}
gated_interaction <- full_validate_distance(gated_interaction, "Gated interaction")
numeric_mean <- full_mean_distinct(numeric_distance)
gamma_values <- c(0, 0.5, 1)
gaussian <- design$spectral$affinities$gaussian
cluster_seed <- artifact$data_seed +
  as.integer(design$execution$cluster_seed_offsets$gaussian)
rows <- lapply(gamma_values, function(gamma) {
  weight <- gamma * numeric_mean / length(categorical_names)
  distance <- full_validate_distance(base + weight * gated_interaction)
  if (gamma == 0 && !identical(distance, base)) {
    stop("Gamma zero does not recover the package baseline exactly.", call. = FALSE)
  }
  scored <- full_score_distance(
    distance, artifact$truth, "gaussian", gaussian, cluster_seed, design
  )
  data.frame(
    index_rule = artifact$index_rule, replicate = artifact$replica,
    n = nrow(data), gamma = gamma, ARI = scored$ARI,
    interaction_weight = weight,
    added_mean = full_mean_distinct(weight * gated_interaction),
    numeric_mean = numeric_mean,
    categorical_mean = full_mean_distinct(categorical_distance),
    gate_seconds = gate_seconds,
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, rows)

sidecar <- list(
  analysis_id = "sim2-exploratory-scale-aware-gated-check",
  note = paste(
    "Paired check on reconstructed, not archived, Simulation 2 data.",
    "Gamma 0 is pure manydist custom pc_scores/TVD/commensurable FALSE.",
    "Gamma 0.5 and 1 are exploratory scale-aware candidates, not locked methods."
  ),
  input_artifact = artifact_path,
  input_artifact_md5 = unname(tools::md5sum(artifact_path)),
  scorer_md5 = unname(tools::md5sum(script)),
  design_config_md5 = attr(design, "config_md5"),
  manydist_source_commit = unname(suppressWarnings(system2(
    "git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE
  ))),
  manydist_version = as.character(utils::packageVersion("manydist")),
  parameters = c(settings, list(prop_nn = prop_nn, score = score,
    decision = decision, gate_permutations = gate$permutations,
    gate_alpha = gate$alpha, gamma_values = gamma_values)),
  decomposition_error = decomposition_error,
  gate_detected = gate$detected,
  gate_p_value = gate$p_value,
  gate_statistic = gate$statistic,
  results = results,
  R_version = R.version.string
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_file <- tempfile("scale-aware-gate-", tmpdir = dirname(output_path),
                           fileext = ".rds")
saveRDS(sidecar, temporary_file, version = 3L)
if (file.exists(output_path)) {
  unlink(temporary_file)
  stop("Sidecar appeared during evaluation; refusing to overwrite it.", call. = FALSE)
}
if (!file.rename(temporary_file, output_path)) {
  unlink(temporary_file)
  stop("Could not atomically save sidecar.", call. = FALSE)
}
cat("Exploratory sidecar completed: ", output_path, "\n", sep = "")
print(results[c("index_rule", "replicate", "gamma", "ARI", "added_mean")],
      row.names = FALSE)
print(gate$detected)
