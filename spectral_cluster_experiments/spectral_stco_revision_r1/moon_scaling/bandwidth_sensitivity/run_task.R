#!/usr/bin/env Rscript

task_script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
}

task_parse_arguments <- function(arguments) {
  output <- list(configuration = NULL, replicate = NULL)
  position <- 1L
  while (position <= length(arguments)) {
    if (position == length(arguments)) stop("Missing argument value.", call. = FALSE)
    key <- arguments[[position]]
    value <- arguments[[position + 1L]]
    if (!key %in% c("--configuration", "--replicate")) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    output[[sub("^--", "", key)]] <- value
    position <- position + 2L
  }
  if (any(vapply(output, is.null, logical(1)))) {
    stop("Use --configuration ID --replicate ID.", call. = FALSE)
  }
  output$replicate <- suppressWarnings(as.integer(output$replicate))
  if (!is.finite(output$replicate) || output$replicate < 1L) {
    stop("`--replicate` must be a positive integer.", call. = FALSE)
  }
  output
}

script <- task_script_path()
bandwidth_dir <- dirname(script)
suite_dir <- dirname(bandwidth_dir)
revision_dir <- dirname(suite_dir)
repo_root <- dirname(dirname(revision_dir))
full_dir <- file.path(revision_dir, "full")
arguments <- task_parse_arguments(commandArgs(trailingOnly = TRUE))

for (package in c("devtools", "digest", "mclust", "yaml")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The bandwidth task requires `", package, "`.", call. = FALSE)
  }
}
devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)

source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))
source(file.path(full_dir, "R", "engine.R"))
source(file.path(suite_dir, "R", "generator.R"))
source(file.path(suite_dir, "R", "gated_distance.R"))
source(file.path(bandwidth_dir, "R", "gaussian_affinity.R"))

design_path <- file.path(bandwidth_dir, "design.yml")
parent_design_path <- file.path(suite_dir, "design.yml")
design <- yaml::read_yaml(design_path)
parent_design <- yaml::read_yaml(parent_design_path)
parent_design_md5 <- unname(tools::md5sum(parent_design_path))
if (!identical(parent_design_md5, as.character(design$parent$design_md5)) ||
    !identical(
      as.character(parent_design$analysis_id),
      as.character(design$parent$analysis_id)
    )) {
  stop("The parent moon-scaling design does not match the sidecar contract.",
       call. = FALSE)
}

configurations <- parent_design$generator$configurations
if (!arguments$configuration %in% names(configurations)) {
  stop("Unknown configuration: ", arguments$configuration, call. = FALSE)
}
if (arguments$replicate > as.integer(design$execution$smoke_replicates)) {
  stop("This runner is restricted to the declared smoke replicates.", call. = FALSE)
}

expected_version <- as.character(parent_design$manydist$version)
installed_version <- as.character(utils::packageVersion("manydist"))
if (!identical(installed_version, expected_version)) {
  stop("Expected manydist ", expected_version, " but found ", installed_version, ".",
       call. = FALSE)
}
source_commit <- suppressWarnings(system2(
  "git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE
))
source_commit <- if (length(source_commit)) source_commit[[1L]] else NA_character_
if (!identical(source_commit, as.character(parent_design$manydist$source_commit))) {
  stop("The manydist source commit does not match the parent manifest.",
       call. = FALSE)
}
if (isTRUE(design$execution$require_clean_manydist_R)) {
  package_changes <- suppressWarnings(system2(
    "git",
    c("-C", repo_root, "status", "--porcelain", "--", "manydist/R"),
    stdout = TRUE,
    stderr = TRUE
  ))
  if (length(package_changes)) {
    stop("The bandwidth sidecar requires a clean `manydist/R` tree.",
         call. = FALSE)
  }
}

manydist_r_files <- sort(list.files(
  file.path(repo_root, "manydist", "R"),
  pattern = "[.]R$", full.names = TRUE
))
source_files <- c(
  design = design_path,
  runner = script,
  affinity = file.path(bandwidth_dir, "R", "gaussian_affinity.R"),
  parent_design = parent_design_path,
  generator = file.path(suite_dir, "R", "generator.R"),
  gated_distance = file.path(suite_dir, "R", "gated_distance.R"),
  bridge = file.path(revision_dir, "R", "manydist_bridge.R"),
  gate_helpers = file.path(revision_dir, "R", "distance_r1.R"),
  spectral_engine = file.path(full_dir, "R", "engine.R"),
  stats::setNames(manydist_r_files, paste0("manydist_R/", basename(manydist_r_files)))
)
source_hashes <- unname(tools::md5sum(source_files))
names(source_hashes) <- names(source_files)
config_tag <- substr(digest::digest(source_hashes, algo = "md5"), 1L, 10L)
output_dir <- file.path(
  bandwidth_dir, "results", "smoke", config_tag, arguments$configuration
)
output_file <- file.path(
  output_dir, sprintf("replicate-%02d.rds", arguments$replicate)
)
if (file.exists(output_file)) {
  cat("Bandwidth artifact already exists; leaving it unchanged: ",
      output_file, "\n", sep = "")
  quit(status = 0L)
}

parent_tag <- substr(parent_design_md5, 1L, 10L)
parent_artifact_path <- file.path(
  suite_dir, "results", "smoke", parent_tag,
  arguments$configuration, as.character(design$parent$interaction_state),
  sprintf("replicate-%02d.rds", arguments$replicate)
)
if (!file.exists(parent_artifact_path)) {
  stop("The required parent smoke artifact is missing: ", parent_artifact_path,
       call. = FALSE)
}
parent_artifact <- readRDS(parent_artifact_path)
if (!identical(parent_artifact$analysis_id, as.character(design$parent$analysis_id)) ||
    !identical(parent_artifact$design_tag, parent_tag)) {
  stop("The parent smoke artifact does not match its declared design.",
       call. = FALSE)
}

configuration <- configurations[[arguments$configuration]]
data_seed <- as.integer(parent_design$execution$data_seed_base) +
  (arguments$replicate - 1L) *
  as.integer(parent_design$execution$replicate_seed_stride)
generated <- moon_scaling_generate(
  seed = data_seed,
  n_per_moon = as.integer(parent_design$generator$n_per_moon),
  p_numeric = as.integer(configuration$p_numeric),
  q_categorical = as.integer(configuration$q_categorical),
  interaction_state = as.character(design$parent$interaction_state),
  moon_noise_sd = as.numeric(parent_design$generator$moon_noise_sd),
  sensor_noise_sd = as.numeric(parent_design$generator$sensor_noise_sd),
  band_probabilities = unlist(parent_design$generator$band_probabilities),
  p_numeric_max = as.integer(parent_design$generator$p_numeric_max),
  q_categorical_max = as.integer(parent_design$generator$q_categorical_max)
)

gated_specification <- parent_design$distance$gated_extension
gate_specification <- gated_specification$gate
gated_started <- proc.time()[[3L]]
gated_fit <- moon_scaling_gated_distance(
  data = generated$data,
  gamma = as.numeric(gated_specification$gamma),
  prop_nn = as.numeric(gated_specification$prop_nn),
  score = as.character(gated_specification$score),
  decision = as.character(gated_specification$decision),
  gate_permutations = as.integer(gate_specification$permutations),
  gate_alpha = as.numeric(gate_specification$alpha),
  gate_seed = data_seed + as.integer(gate_specification$seed_offset)
)
gated_seconds <- proc.time()[[3L]] - gated_started

distance_by_method <- list(
  ab_gated = gated_fit$distance,
  ab_no_interaction = gated_fit$distance_no_interaction
)
distance_seconds <- c(
  ab_gated = gated_seconds,
  ab_no_interaction = gated_fit$diagnostics$baseline_seconds
)
competitor_presets <- vapply(
  parent_design$distance$competitors,
  function(x) as.character(x$preset),
  character(1)
)
for (method in names(competitor_presets)) {
  started <- proc.time()[[3L]]
  fit <- manydist::mdist(generated$data, preset = competitor_presets[[method]])
  distance_by_method[[method]] <- full_validate_distance(
    as.matrix(fit$distance), method
  )
  distance_seconds[[method]] <- proc.time()[[3L]] - started
}

multipliers <- as.numeric(unlist(design$bandwidth$multipliers))
primary_multiplier <- as.numeric(design$bandwidth$primary_multiplier)
cluster_seed <- data_seed +
  as.integer(parent_design$execution$cluster_seed_offsets$gaussian)
kmeans_specification <- parent_design$spectral$kmeans
result_rows <- list()
row_index <- 1L
for (method in names(distance_by_method)) {
  for (multiplier in multipliers) {
    started <- proc.time()[[3L]]
    score <- moon_bandwidth_score_distance(
      distance = distance_by_method[[method]],
      truth = generated$truth,
      multiplier = multiplier,
      cluster_seed = cluster_seed,
      kmeans_specification = kmeans_specification
    )
    result_rows[[row_index]] <- data.frame(
      configuration = arguments$configuration,
      interaction_state = as.character(design$parent$interaction_state),
      replicate = arguments$replicate,
      data_seed = data_seed,
      n = nrow(generated$data),
      p_numeric = as.integer(configuration$p_numeric),
      q_categorical = as.integer(configuration$q_categorical),
      method = method,
      affinity = "gaussian",
      bandwidth_base_rule = as.character(design$bandwidth$base_rule),
      bandwidth_multiplier = multiplier,
      bandwidth_role = if (multiplier == primary_multiplier) "primary" else "sensitivity",
      affinity_base_scale = score$affinity_base_scale,
      affinity_scale = score$affinity_scale,
      ARI = score$ARI,
      ifault = score$ifault,
      distance_seconds = unname(distance_seconds[[method]]),
      clustering_seconds = proc.time()[[3L]] - started,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
results <- do.call(rbind, result_rows)

parent_primary <- parent_artifact$results[
  parent_artifact$results$affinity == "gaussian",
  c("method", "ARI")
]
current_primary <- results[
  results$bandwidth_multiplier == primary_multiplier,
  c("method", "ARI")
]
primary_check <- merge(
  parent_primary, current_primary,
  by = "method", suffixes = c("_parent", "_sidecar"), sort = TRUE
)
if (nrow(primary_check) != length(distance_by_method) ||
    any(abs(primary_check$ARI_parent - primary_check$ARI_sidecar) > 1e-12)) {
  stop("The sidecar does not reproduce the parent Gaussian results at c = 1.",
       call. = FALSE)
}

artifact <- list(
  analysis_id = as.character(design$analysis_id),
  status = as.character(design$status),
  parent_analysis_id = as.character(design$parent$analysis_id),
  parent_artifact = parent_artifact_path,
  parent_primary_check = primary_check,
  configuration = arguments$configuration,
  interaction_state = as.character(design$parent$interaction_state),
  replicate = arguments$replicate,
  data_seed = data_seed,
  design_tag = config_tag,
  source_hashes = source_hashes,
  source_commit = source_commit,
  manydist_version = installed_version,
  results = results,
  gate = parent_artifact$gate,
  R_version = R.version.string,
  package_versions = c(
    manydist = installed_version,
    digest = as.character(utils::packageVersion("digest")),
    mclust = as.character(utils::packageVersion("mclust")),
    yaml = as.character(utils::packageVersion("yaml"))
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
temporary_file <- tempfile("sim3-bandwidth-", tmpdir = output_dir, fileext = ".rds")
saveRDS(artifact, temporary_file, version = 3L)
if (file.exists(output_file)) {
  unlink(temporary_file)
  stop("Output appeared during the task; refusing to overwrite it.", call. = FALSE)
}
if (!file.rename(temporary_file, output_file)) {
  unlink(temporary_file)
  stop("Could not atomically save the bandwidth artifact.", call. = FALSE)
}

cat("Bandwidth smoke task completed: ", output_file, "\n", sep = "")
print(results[c(
  "method", "bandwidth_multiplier", "ARI", "affinity_scale"
)], row.names = FALSE)

