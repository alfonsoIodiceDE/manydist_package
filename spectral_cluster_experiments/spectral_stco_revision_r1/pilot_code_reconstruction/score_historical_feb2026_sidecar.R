#!/usr/bin/env Rscript

# Score the exact February-2026 source formula on an existing paired-pilot
# dataset. Output lives separately from the immutable original pilot artifact.
# This uses the *current pilot's* common spectral engine/affinities so results
# compare fairly with its R1 methods; it does not recreate the old full run.

argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
pilot_dir <- dirname(script)
revision_dir <- dirname(pilot_dir)
full_dir <- file.path(revision_dir, "full")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript score_historical_feb2026_sidecar.R PATH_TO_PILOT_RDS",
       call. = FALSE)
}
artifact_path <- normalizePath(arguments[[1L]], mustWork = TRUE)
results_dir <- paste0(normalizePath(file.path(pilot_dir, "results"), mustWork = TRUE),
                      .Platform$file.sep)
if (!startsWith(artifact_path, results_dir) ||
    !grepl("\\.rds$", artifact_path)) {
  stop("The input must be an existing pilot RDS under results/.", call. = FALSE)
}
relative_path <- substring(artifact_path, nchar(results_dir) + 1L)
# Version 1 sidecars used a misoriented categorical profile normalization.
# Preserve them for auditability and write corrected results to a new tree.
output_path <- file.path(pilot_dir, "historical_feb2026_sidecars_v2", relative_path)
if (file.exists(output_path)) {
  cat("Historical sidecar already exists; leaving it unchanged: ", output_path,
      "\n", sep = "")
  quit(status = 0L)
}

source(file.path(full_dir, "R", "design.R"))
source(file.path(full_dir, "R", "engine.R"))
source(file.path(pilot_dir, "R", "historical_feb2026_distance.R"))
if (!requireNamespace("mclust", quietly = TRUE)) {
  stop("mclust is required for scoring.", call. = FALSE)
}
design <- full_read_design(file.path(full_dir, "design.yml"))
artifact <- readRDS(artifact_path)
if (!identical(artifact$analysis_id, "sim2-code-reconstruction-paired-pilot") ||
    !identical(artifact$design_config_md5, attr(design, "config_md5"))) {
  stop("Input is not a paired pilot artifact using this design.", call. = FALSE)
}

start <- proc.time()[[3L]]
distance <- historical_feb2026_sim2_distance(artifact$data)
distance_seconds <- proc.time()[[3L]] - start
affinities <- full_enabled_affinities(design)
rows <- lapply(names(affinities), function(affinity_name) {
  cluster_seed <- as.integer(artifact$data_seed) +
    as.integer(design$execution$cluster_seed_offsets[[affinity_name]])
  cluster_start <- proc.time()[[3L]]
  score <- full_score_distance(
    distance, artifact$truth, affinity_name, affinities[[affinity_name]],
    cluster_seed, design
  )
  data.frame(
    replicate = artifact$replica,
    data_seed = artifact$data_seed,
    index_rule = artifact$index_rule,
    n = nrow(artifact$data),
    method = "historical_feb2026_euclidean_pc_tvd_noop_v2",
    affinity = affinity_name,
    ARI = score$ARI,
    ifault = score$ifault,
    distance_seconds = distance_seconds,
    clustering_seconds = proc.time()[[3L]] - cluster_start,
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, rows)
sidecar <- list(
  analysis_id = "sim2-code-reconstruction-historical-feb2026-sidecar-v2",
  comparator_version = 2L,
  historical_source_commit = "79563a7",
  historical_source_branch = paste(
    "mdist custom: Euclidean numeric PC scores + categorical total variation;",
    "interaction ignored in this branch"
  ),
  caveat = paste(
    "February mdist() lacked prop_nn and score arguments. The saved April",
    "simulation xBIG.R call cannot run verbatim against this version; exact",
    "archived-run package version/driver remain unverified. This sidecar uses",
    "the current paired pilot's common spectral engine, not old clustering."
  ),
  input_artifact = artifact_path,
  input_artifact_md5 = unname(tools::md5sum(artifact_path)),
  comparator_source_md5 = unname(tools::md5sum(file.path(
    pilot_dir, "R", "historical_feb2026_distance.R"
  ))),
  design_config_md5 = attr(design, "config_md5"),
  results = results,
  mean_distance = full_mean_distinct(distance),
  R_version = R.version.string
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_file <- tempfile("historical-", tmpdir = dirname(output_path),
                           fileext = ".rds")
saveRDS(sidecar, temporary_file, version = 3L)
if (file.exists(output_path)) {
  unlink(temporary_file)
  stop("Sidecar appeared during evaluation; refusing to overwrite it.",
       call. = FALSE)
}
if (!file.rename(temporary_file, output_path)) {
  unlink(temporary_file)
  stop("Could not atomically save sidecar.", call. = FALSE)
}
cat("Historical sidecar completed: ", output_path, "\n", sep = "")
print(results[c("method", "affinity", "ARI", "distance_seconds")], row.names = FALSE)
