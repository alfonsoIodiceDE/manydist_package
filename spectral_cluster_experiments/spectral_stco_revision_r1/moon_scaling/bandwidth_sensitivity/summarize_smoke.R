#!/usr/bin/env Rscript

script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
}

script <- script_path()
bandwidth_dir <- dirname(script)
suite_dir <- dirname(bandwidth_dir)
revision_dir <- dirname(suite_dir)
repo_root <- dirname(dirname(revision_dir))
full_dir <- file.path(revision_dir, "full")
for (package in c("digest", "yaml")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The bandwidth summarizer requires `", package, "`.", call. = FALSE)
  }
}

design_path <- file.path(bandwidth_dir, "design.yml")
parent_design_path <- file.path(suite_dir, "design.yml")
design <- yaml::read_yaml(design_path)
parent_design <- yaml::read_yaml(parent_design_path)
manydist_r_files <- sort(list.files(
  file.path(repo_root, "manydist", "R"),
  pattern = "[.]R$", full.names = TRUE
))
source_files <- c(
  design = design_path,
  runner = file.path(bandwidth_dir, "run_task.R"),
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
design_tag <- substr(digest::digest(source_hashes, algo = "md5"), 1L, 10L)
input_dir <- file.path(bandwidth_dir, "results", "smoke", design_tag)

configurations <- names(parent_design$generator$configurations)
replicates <- seq_len(as.integer(design$execution$smoke_replicates))
expected <- expand.grid(
  configuration = configurations,
  replicate = replicates,
  stringsAsFactors = FALSE
)
expected$file <- file.path(
  input_dir, expected$configuration,
  sprintf("replicate-%02d.rds", expected$replicate)
)
missing <- expected$file[!file.exists(expected$file)]
if (length(missing)) {
  stop("Missing expected artifact(s):\n", paste(missing, collapse = "\n"),
       call. = FALSE)
}

artifacts <- lapply(expected$file, readRDS)
valid <- vapply(seq_along(artifacts), function(index) {
  artifact <- artifacts[[index]]
  identical(artifact$analysis_id, as.character(design$analysis_id)) &&
    identical(artifact$design_tag, design_tag) &&
    identical(artifact$configuration, expected$configuration[[index]]) &&
    identical(as.integer(artifact$replicate), expected$replicate[[index]]) &&
    identical(artifact$source_hashes, source_hashes)
}, logical(1))
if (!all(valid)) {
  stop("At least one artifact does not match the current frozen sources.",
       call. = FALSE)
}

metrics <- do.call(rbind, lapply(artifacts, `[[`, "results"))
rownames(metrics) <- NULL
methods <- c(
  "ab_gated", "ab_no_interaction", "gower", "modified_gower",
  "standardized_euclidean_onehot"
)
multipliers <- as.numeric(unlist(design$bandwidth$multipliers))
metrics <- metrics[order(
  match(metrics$configuration, configurations),
  match(metrics$bandwidth_multiplier, multipliers),
  match(metrics$method, methods)
), ]

comparisons <- do.call(rbind, lapply(
  split(seq_len(nrow(metrics)), interaction(
    metrics$configuration, metrics$bandwidth_multiplier, drop = TRUE
  )),
  function(rows) {
    block <- metrics[rows, ]
    gated <- block$ARI[block$method == "ab_gated"]
    baseline <- block$ARI[block$method == "ab_no_interaction"]
    competitors <- block[!block$method %in% c("ab_gated", "ab_no_interaction"), ]
    best_row <- competitors[which.max(competitors$ARI), ]
    data.frame(
      configuration = block$configuration[[1L]],
      bandwidth_multiplier = block$bandwidth_multiplier[[1L]],
      ARI_gated = gated,
      ARI_no_interaction = baseline,
      gated_minus_no_interaction = gated - baseline,
      best_external_method = best_row$method,
      best_external_ARI = best_row$ARI,
      gated_minus_best_external = gated - best_row$ARI,
      gated_is_best = gated >= max(block$ARI),
      stringsAsFactors = FALSE
    )
  }
))
comparisons <- comparisons[order(
  match(comparisons$configuration, configurations),
  match(comparisons$bandwidth_multiplier, multipliers)
), ]
rownames(comparisons) <- NULL

write.csv(metrics, file.path(input_dir, "metrics.csv"), row.names = FALSE)
write.csv(comparisons, file.path(input_dir, "comparisons.csv"), row.names = FALSE)
cat("Validated and summarized ", length(artifacts), " artifact(s).\n", sep = "")
cat("Design tag: ", design_tag, "\n", sep = "")
cat("Output: ", input_dir, "\n", sep = "")
print(comparisons, row.names = FALSE)

