#!/usr/bin/env Rscript

# Paired original-size pilot on a reconstruction of the available historical
# Simulation 2 script. This does not alter the frozen full-study manifest.

pilot_script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
}

pilot_parse_args <- function(arguments) {
  out <- list(replicate = NULL, index_rule = "literal", n = 500L)
  position <- 1L
  while (position <= length(arguments)) {
    key <- arguments[[position]]
    if (!key %in% c("--replicate", "--index-rule", "--n") ||
        position == length(arguments)) {
      stop("Expected --replicate ID [--index-rule literal|corrected] [--n SIZE].",
           call. = FALSE)
    }
    value <- arguments[[position + 1L]]
    out[[gsub("-", "_", sub("^--", "", key))]] <- value
    position <- position + 2L
  }
  out$replicate <- suppressWarnings(as.integer(out$replicate))
  out$n <- suppressWarnings(as.integer(out$n))
  if (length(out$replicate) != 1L || !is.finite(out$replicate) ||
      out$replicate < 1L || out$replicate > 50L) {
    stop("--replicate must be an integer from 1 to 50.", call. = FALSE)
  }
  if (length(out$n) != 1L || !is.finite(out$n) ||
      out$n < 100L || out$n %% 10L != 0L) {
    stop("--n must be an integer >= 100 divisible by 10.", call. = FALSE)
  }
  if (!out$index_rule %in% c("literal", "corrected")) {
    stop("--index-rule must be literal or corrected.", call. = FALSE)
  }
  out
}

script <- pilot_script_path()
pilot_dir <- dirname(script)
revision_dir <- dirname(pilot_dir)
repo_root <- dirname(dirname(revision_dir))
full_dir <- file.path(revision_dir, "full")
arguments <- pilot_parse_args(commandArgs(trailingOnly = TRUE))

source(file.path(full_dir, "R", "design.R"))
source(file.path(pilot_dir, "R", "sim2_generator.R"))
source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))
source(file.path(full_dir, "R", "engine.R"))

design <- full_read_design(file.path(full_dir, "design.yml"))
config_tag <- substr(attr(design, "config_md5"), 1L, 10L)
output_dir <- file.path(
  pilot_dir, "results", config_tag,
  sprintf("n%04d", arguments$n), arguments$index_rule
)
output_file <- file.path(
  output_dir, sprintf("replicate-%02d.rds", arguments$replicate)
)
if (file.exists(output_file)) {
  cat("Pilot artifact already exists; leaving it unchanged: ", output_file, "\n", sep = "")
  quit(status = 0L)
}

for (package in c("devtools", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The pilot requires the package ", package, ".", call. = FALSE)
  }
}
devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)

data_seed <- as.integer(26091900L + 100L * arguments$replicate)
generated <- sim2_reconstructed_data(
  data_seed,
  n = arguments$n,
  index_rule = arguments$index_rule
)
started <- Sys.time()

r1_start <- proc.time()[[3L]]
r1 <- scmix_distance_r1(
  data = generated$x,
  prop_nn = as.numeric(design$distance$defaults$prop_nn),
  score = as.character(design$distance$defaults$score),
  decision = as.character(design$distance$defaults$decision),
  interaction_gate = as.character(design$distance$categorical$interaction_gate$method),
  gate_permutations = as.integer(design$distance$categorical$interaction_gate$permutations),
  gate_alpha = as.numeric(design$distance$categorical$interaction_gate$alpha),
  gate_seed = data_seed +
    as.integer(design$distance$categorical$interaction_gate$seed_offset),
  strict_provenance = isTRUE(design$manydist$strict_provenance),
  manydist_repo = repo_root
)
r1_seconds <- proc.time()[[3L]] - r1_start

# Closest available current-package diagnostic to the submitted Sim2 call.
# It retains the historical prop_nn and score, but current manydist custom
# numerical geometry is not guaranteed identical to the old Euclidean route.
legacy_start <- proc.time()[[3L]]
legacy <- manydist::mdist(
  generated$x,
  preset = "custom",
  method_cat = "tvd",
  method_num = "pc_scores",
  commensurable = FALSE,
  interaction = TRUE,
  prop_nn = 0.05,
  score = "logloss"
)
legacy_seconds <- proc.time()[[3L]] - legacy_start

distances <- list(
  AB_BW_R1_interaction = as.matrix(r1$distance),
  AB_BW_R1_no_interaction = as.matrix(r1$distance_no_interaction),
  legacy_current_custom_logloss = as.matrix(legacy$distance)
)
distance_seconds <- c(
  AB_BW_R1_interaction = r1_seconds,
  AB_BW_R1_no_interaction = r1_seconds,
  legacy_current_custom_logloss = legacy_seconds
)
affinities <- full_enabled_affinities(design)

results <- do.call(rbind, lapply(names(distances), function(method) {
  do.call(rbind, lapply(names(affinities), function(affinity) {
    cluster_seed <- data_seed +
      as.integer(design$execution$cluster_seed_offsets[[affinity]])
    cluster_start <- proc.time()[[3L]]
    score <- full_score_distance(
      distances[[method]],
      generated$truth,
      affinity,
      affinities[[affinity]],
      cluster_seed,
      design
    )
    data.frame(
      replicate = arguments$replicate,
      data_seed = data_seed,
      index_rule = arguments$index_rule,
      n = arguments$n,
      method = method,
      affinity = affinity,
      ARI = score$ARI,
      ifault = score$ifault,
      distance_seconds = unname(distance_seconds[[method]]),
      clustering_seconds = proc.time()[[3L]] - cluster_start,
      stringsAsFactors = FALSE
    )
  }))
}))

gate <- r1$interaction_gate
gate_diagnostics <- data.frame(
  variable = names(gate$detected),
  detected = unname(gate$detected),
  p_value = unname(gate$p_value),
  statistic = unname(gate$statistic),
  stringsAsFactors = FALSE
)
component_diagnostics <- list(
  numeric_mean = full_mean_distinct(r1$components$numeric$scaled),
  categorical_tvd_mean = full_mean_distinct(r1$components$categorical_tvd_sum),
  interaction_gated_mean = full_mean_distinct(
    r1$components$categorical_interaction_gated_sum
  ),
  rho = r1$config$rho
)

artifact <- list(
  analysis_id = "sim2-code-reconstruction-paired-pilot",
  source_script = script,
  source_script_md5 = unname(tools::md5sum(script)),
  generator_script_md5 = unname(tools::md5sum(
    file.path(pilot_dir, "R", "sim2_generator.R")
  )),
  design_config_md5 = attr(design, "config_md5"),
  replica = arguments$replicate,
  index_rule = arguments$index_rule,
  data_seed = data_seed,
  data = generated$x,
  truth = generated$truth,
  generator_diagnostics = generated$diagnostics,
  results = results,
  gate_diagnostics = gate_diagnostics,
  component_diagnostics = component_diagnostics,
  r1_provenance = r1$provenance,
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  R_version = R.version.string,
  package_versions = c(
    manydist = as.character(utils::packageVersion("manydist")),
    mclust = as.character(utils::packageVersion("mclust"))
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
temporary_file <- tempfile("replicate-", tmpdir = output_dir, fileext = ".rds")
saveRDS(artifact, temporary_file, version = 3L)
if (file.exists(output_file)) {
  unlink(temporary_file)
  stop("Pilot output appeared during evaluation; refusing to overwrite it.", call. = FALSE)
}
if (!file.rename(temporary_file, output_file)) {
  unlink(temporary_file)
  stop("Could not atomically save the pilot artifact.", call. = FALSE)
}

cat("Reconstruction pilot completed: ", output_file, "\n", sep = "")
cat("Elapsed seconds: ", sprintf("%.2f", artifact$elapsed_seconds), "\n", sep = "")
print(results[c("method", "affinity", "ARI", "distance_seconds")], row.names = FALSE)
