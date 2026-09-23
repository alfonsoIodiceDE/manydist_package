#!/usr/bin/env Rscript

# Evaluate the author's requested current-manydist custom interaction switch
# on an existing paired Simulation 2 pilot dataset. This is a package-native
# configuration comparison, not an interaction-only ablation: current mdist()
# also changes the categorical block multiplier when interaction is TRUE.

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE)
pilot_dir <- dirname(script)
revision_dir <- dirname(pilot_dir)
repo_root <- dirname(dirname(revision_dir))
full_dir <- file.path(revision_dir, "full")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript score_custom_mdist_pair.R PATH_TO_PILOT_RDS", call. = FALSE)
}

artifact_path <- normalizePath(arguments[[1L]], mustWork = TRUE)
results_root <- paste0(normalizePath(file.path(pilot_dir, "results"), mustWork = TRUE),
                       .Platform$file.sep)
if (!startsWith(artifact_path, results_root) ||
    !grepl("[.]rds$", artifact_path)) {
  stop("Input must be a saved pilot RDS under results/.", call. = FALSE)
}
relative_path <- substring(artifact_path, nchar(results_root) + 1L)
output_path <- file.path(pilot_dir, "custom_mdist_pair_sidecars", relative_path)
if (file.exists(output_path)) {
  cat("Custom mdist sidecar already exists; leaving it unchanged: ",
      output_path, "\n", sep = "")
  quit(status = 0L)
}

source(file.path(full_dir, "R", "design.R"))
source(file.path(full_dir, "R", "engine.R"))
if (!requireNamespace("devtools", quietly = TRUE) ||
    !requireNamespace("mclust", quietly = TRUE)) {
  stop("This scorer requires devtools and mclust.", call. = FALSE)
}
devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)

design <- full_read_design(file.path(full_dir, "design.yml"))
artifact <- readRDS(artifact_path)
if (!identical(artifact$analysis_id, "sim2-code-reconstruction-paired-pilot") ||
    !identical(artifact$design_config_md5, attr(design, "config_md5"))) {
  stop("Input is not a paired pilot artifact using this design.", call. = FALSE)
}

prop_nn <- as.numeric(design$distance$defaults$prop_nn)
score <- as.character(design$distance$defaults$score)
decision <- as.character(design$distance$defaults$decision)
affinities <- full_enabled_affinities(design)
distance_by_setting <- list()
distance_seconds <- numeric()

for (setting in c("off", "on")) {
  start <- proc.time()[[3L]]
  fit <- manydist::mdist(
    artifact$data,
    preset = "custom",
    method_num = "pc_scores",
    method_cat = "tvd",
    commensurable = FALSE,
    interaction = identical(setting, "on"),
    prop_nn = prop_nn,
    score = score,
    decision = decision
  )
  distance_by_setting[[setting]] <- full_validate_distance(
    as.matrix(fit$distance), paste("Custom mdist", setting)
  )
  distance_seconds[[setting]] <- proc.time()[[3L]] - start
}

rows <- lapply(names(distance_by_setting), function(setting) {
  do.call(rbind, lapply(names(affinities), function(affinity) {
    cluster_seed <- artifact$data_seed +
      as.integer(design$execution$cluster_seed_offsets[[affinity]])
    start <- proc.time()[[3L]]
    fit <- full_score_distance(
      distance_by_setting[[setting]], artifact$truth, affinity,
      affinities[[affinity]], cluster_seed, design
    )
    data.frame(
      replicate = artifact$replica,
      index_rule = artifact$index_rule,
      n = nrow(artifact$data),
      method = paste0("mdist_custom_pc_tvd_comm_false_interaction_", setting),
      affinity = affinity,
      ARI = fit$ARI,
      ifault = fit$ifault,
      distance_seconds = unname(distance_seconds[[setting]]),
      clustering_seconds = proc.time()[[3L]] - start,
      stringsAsFactors = FALSE
    )
  }))
})
results <- do.call(rbind, rows)

off <- distance_by_setting$off
on <- distance_by_setting$on
sidecar <- list(
  analysis_id = "sim2-code-reconstruction-current-mdist-custom-pair",
  note = paste(
    "This is the author's requested package-native comparison. In the",
    "current mdist implementation, interaction TRUE also multiplies the",
    "categorical block by q_cat/q_num, so it is not an interaction-only",
    "ablation or the earlier gated u_dep_bw method."
  ),
  input_artifact = artifact_path,
  input_artifact_md5 = unname(tools::md5sum(artifact_path)),
  scorer_md5 = unname(tools::md5sum(script)),
  design_config_md5 = attr(design, "config_md5"),
  manydist_source_commit = unname(suppressWarnings(system2(
    "git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE
  ))),
  manydist_version = as.character(utils::packageVersion("manydist")),
  parameters = list(
    preset = "custom", method_num = "pc_scores", method_cat = "tvd",
    commensurable = FALSE, prop_nn = prop_nn, score = score,
    decision = decision
  ),
  mean_distance = c(off = full_mean_distinct(off),
                    on = full_mean_distinct(on)),
  mean_abs_distance_change = mean(abs(on[upper.tri(on)] - off[upper.tri(off)])),
  max_abs_distance_change = max(abs(on - off)),
  results = results,
  R_version = R.version.string
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_file <- tempfile("custom-mdist-", tmpdir = dirname(output_path),
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
cat("Custom mdist pair completed: ", output_path, "\n", sep = "")
print(results[c("method", "affinity", "ARI", "distance_seconds")],
      row.names = FALSE)
