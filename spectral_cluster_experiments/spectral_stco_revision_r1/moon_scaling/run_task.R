#!/usr/bin/env Rscript

task_script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
}

task_parse_arguments <- function(arguments) {
  output <- list(configuration = NULL, state = NULL, replicate = NULL)
  position <- 1L
  while (position <= length(arguments)) {
    if (position == length(arguments)) stop("Missing argument value.", call. = FALSE)
    key <- arguments[[position]]
    value <- arguments[[position + 1L]]
    if (!key %in% c("--configuration", "--state", "--replicate")) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    output[[sub("^--", "", key)]] <- value
    position <- position + 2L
  }
  if (any(vapply(output, is.null, logical(1)))) {
    stop("Use --configuration ID --state signal|null --replicate ID.", call. = FALSE)
  }
  output$replicate <- suppressWarnings(as.integer(output$replicate))
  if (!is.finite(output$replicate) || output$replicate < 1L) {
    stop("`--replicate` must be a positive integer.", call. = FALSE)
  }
  output
}

script <- task_script_path()
suite_dir <- dirname(script)
revision_dir <- dirname(suite_dir)
repo_root <- dirname(dirname(revision_dir))
full_dir <- file.path(revision_dir, "full")
arguments <- task_parse_arguments(commandArgs(trailingOnly = TRUE))

for (package in c("devtools", "mclust", "yaml")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The Simulation 3 task requires `", package, "`.", call. = FALSE)
  }
}
devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)

source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))
source(file.path(full_dir, "R", "engine.R"))
source(file.path(suite_dir, "R", "generator.R"))
source(file.path(suite_dir, "R", "gated_distance.R"))

design_path <- file.path(suite_dir, "design.yml")
design <- yaml::read_yaml(design_path)
configurations <- design$generator$configurations
if (!arguments$configuration %in% names(configurations)) {
  stop("Unknown configuration: ", arguments$configuration, call. = FALSE)
}
if (!arguments$state %in% unlist(design$generator$interaction_states)) {
  stop("Unknown interaction state: ", arguments$state, call. = FALSE)
}
if (arguments$replicate > as.integer(design$execution$smoke_replicates)) {
  stop("This runner is currently restricted to declared smoke replicates.", call. = FALSE)
}

expected_version <- as.character(design$manydist$version)
installed_version <- as.character(utils::packageVersion("manydist"))
if (!identical(installed_version, expected_version)) {
  stop("Expected manydist ", expected_version, " but found ", installed_version, ".",
       call. = FALSE)
}
source_commit <- suppressWarnings(system2(
  "git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE
))
source_commit <- if (length(source_commit)) source_commit[[1L]] else NA_character_
if (!identical(source_commit, as.character(design$manydist$source_commit))) {
  stop("The manydist source commit does not match the Simulation 3 manifest.",
       call. = FALSE)
}

source_files <- c(
  design = design_path,
  runner = script,
  generator = file.path(suite_dir, "R", "generator.R"),
  gated_distance = file.path(suite_dir, "R", "gated_distance.R"),
  bridge = file.path(revision_dir, "R", "manydist_bridge.R"),
  gate_helpers = file.path(revision_dir, "R", "distance_r1.R"),
  spectral_engine = file.path(full_dir, "R", "engine.R")
)
source_hashes <- unname(tools::md5sum(source_files))
names(source_hashes) <- names(source_files)
config_tag <- substr(source_hashes[["design"]], 1L, 10L)
output_dir <- file.path(
  suite_dir, "results", "smoke", config_tag,
  arguments$configuration, arguments$state
)
output_file <- file.path(
  output_dir, sprintf("replicate-%02d.rds", arguments$replicate)
)
if (file.exists(output_file)) {
  cat("Simulation 3 artifact already exists; leaving it unchanged: ",
      output_file, "\n", sep = "")
  quit(status = 0L)
}

configuration <- configurations[[arguments$configuration]]
data_seed <- as.integer(design$execution$data_seed_base) +
  (arguments$replicate - 1L) * as.integer(design$execution$replicate_seed_stride)
generated <- moon_scaling_generate(
  seed = data_seed,
  n_per_moon = as.integer(design$generator$n_per_moon),
  p_numeric = as.integer(configuration$p_numeric),
  q_categorical = as.integer(configuration$q_categorical),
  interaction_state = arguments$state,
  moon_noise_sd = as.numeric(design$generator$moon_noise_sd),
  sensor_noise_sd = as.numeric(design$generator$sensor_noise_sd),
  band_probabilities = unlist(design$generator$band_probabilities),
  p_numeric_max = as.integer(design$generator$p_numeric_max),
  q_categorical_max = as.integer(design$generator$q_categorical_max)
)

gated_specification <- design$distance$gated_extension
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
  design$distance$competitors,
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

affinities <- design$spectral$affinities
score_design <- list(spectral = list(kmeans = design$spectral$kmeans))
result_rows <- list()
row_index <- 1L
for (method in names(distance_by_method)) {
  for (affinity in names(affinities)) {
    cluster_seed <- data_seed +
      as.integer(design$execution$cluster_seed_offsets[[affinity]])
    started <- proc.time()[[3L]]
    score <- full_score_distance(
      distance = distance_by_method[[method]],
      truth = generated$truth,
      affinity_name = affinity,
      affinity_specification = affinities[[affinity]],
      cluster_seed = cluster_seed,
      design = score_design
    )
    result_rows[[row_index]] <- data.frame(
      configuration = arguments$configuration,
      interaction_state = arguments$state,
      replicate = arguments$replicate,
      data_seed = data_seed,
      n = nrow(generated$data),
      p_numeric = as.integer(configuration$p_numeric),
      q_categorical = as.integer(configuration$q_categorical),
      method = method,
      affinity = affinity,
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

artifact <- list(
  analysis_id = as.character(design$analysis_id),
  status = as.character(design$status),
  configuration = arguments$configuration,
  interaction_state = arguments$state,
  replicate = arguments$replicate,
  data_seed = data_seed,
  design_tag = config_tag,
  source_hashes = source_hashes,
  source_commit = source_commit,
  manydist_version = installed_version,
  generator_diagnostics = generated$diagnostics,
  truth_counts = table(generated$truth),
  results = results,
  gate = list(
    detected = gated_fit$gate$detected,
    p_value = gated_fit$gate$p_value,
    statistic = gated_fit$gate$statistic,
    active_names = gated_fit$diagnostics$active_names,
    active_count = gated_fit$diagnostics$active_count
  ),
  distance_diagnostics = gated_fit$diagnostics,
  configuration_record = gated_fit$configuration,
  R_version = R.version.string,
  package_versions = c(
    manydist = installed_version,
    mclust = as.character(utils::packageVersion("mclust")),
    yaml = as.character(utils::packageVersion("yaml"))
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
temporary_file <- tempfile("sim3-", tmpdir = output_dir, fileext = ".rds")
saveRDS(artifact, temporary_file, version = 3L)
if (file.exists(output_file)) {
  unlink(temporary_file)
  stop("Output appeared during the task; refusing to overwrite it.", call. = FALSE)
}
if (!file.rename(temporary_file, output_file)) {
  unlink(temporary_file)
  stop("Could not atomically save the Simulation 3 artifact.", call. = FALSE)
}

cat("Simulation 3 smoke task completed: ", output_file, "\n", sep = "")
print(results[c("method", "affinity", "ARI", "distance_seconds")], row.names = FALSE)
cat("Gate decisions:\n")
print(gated_fit$gate$detected)
