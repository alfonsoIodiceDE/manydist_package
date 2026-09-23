#!/usr/bin/env Rscript

full_script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(argument)) {
    return(normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE))
  }
  normalizePath(
    "spectral_cluster_experiments/revision_r1/full/merge_results.R",
    mustWork = TRUE
  )
}

script_file <- full_script_path()
full_dir <- dirname(script_file)
source(file.path(full_dir, "R", "design.R"))

arguments <- full_parse_cli(commandArgs(trailingOnly = TRUE), mode = "merge")
config_path <- arguments$config %||% file.path(full_dir, "design.yml")
design <- full_read_design(config_path)
tasks <- full_expand_tasks(design)
smoke_mode <- isTRUE(arguments$include_smoke)
input_dir <- full_resolve_output_dir(design, smoke = smoke_mode)
output_dir <- file.path(
  full_resolve_merged_dir(design),
  if (smoke_mode) "smoke" else "full"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
} else {
  character()
}
if (!length(files)) {
  stop("No task artifacts found in: ", input_dir, call. = FALSE)
}

artifacts <- lapply(files, readRDS)
valid <- vapply(artifacts, function(x) {
  identical(x$analysis_id, design$analysis_id) &&
    identical(x$config_md5, attr(design, "config_md5")) &&
    identical(isTRUE(x$smoke), smoke_mode)
}, logical(1))
if (!all(valid)) {
  stop("At least one task artifact does not match this manifest/hash/mode.", call. = FALSE)
}

task_ids <- vapply(artifacts, function(x) as.integer(x$task$task_id), integer(1))
if (anyDuplicated(task_ids)) {
  stop("Duplicate task identifiers found during merge.", call. = FALSE)
}
ordering <- order(task_ids)
artifacts <- artifacts[ordering]
task_ids <- task_ids[ordering]
files <- files[ordering]

bind_rows <- function(rows) {
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x), logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

all_parameter_names <- sort(unique(unlist(lapply(artifacts, function(x) {
  names(x$effective_parameters %||% list())
}))))
parameter_rows <- lapply(artifacts, function(artifact) {
  values <- artifact$effective_parameters %||% list()
  row <- data.frame(task_id = artifact$task$task_id, stringsAsFactors = FALSE)
  for (name in all_parameter_names) {
    row[[paste0("param_", name)]] <- if (is.null(values[[name]])) {
      NA_character_
    } else {
      full_compact_value(values[[name]])
    }
  }
  row
})
parameters <- bind_rows(parameter_rows)

metrics <- bind_rows(lapply(artifacts, `[[`, "results"))
diagnostics <- bind_rows(lapply(artifacts, `[[`, "interaction_diagnostics"))
if (nrow(metrics)) metrics <- merge(metrics, parameters, by = "task_id", all.x = TRUE, sort = FALSE)
if (nrow(diagnostics)) {
  diagnostics <- merge(diagnostics, parameters, by = "task_id", all.x = TRUE, sort = FALSE)
}

status <- bind_rows(lapply(artifacts, function(x) {
  data.frame(
    task_id = x$task$task_id,
    scenario = x$task$scenario,
    replicate = x$task$replicate,
    status = x$status,
    elapsed_seconds = x$elapsed_seconds,
    warning_count = length(x$warnings),
    warnings = paste(x$warnings, collapse = " | "),
    error = x$error,
    output_file = files[match(x$task$task_id, task_ids)],
    stringsAsFactors = FALSE
  )
}))

method_failures <- if (nrow(metrics)) {
  metrics[!is.na(metrics$error) & nzchar(metrics$error), , drop = FALSE]
} else {
  data.frame()
}

summarise_metrics <- function(metrics) {
  if (!nrow(metrics)) return(data.frame())
  parameter_columns <- grep("^param_", names(metrics), value = TRUE)
  group_columns <- c(
    "scenario", parameter_columns, "setting_id", "prop_nn", "gamma",
    "numeric_block_weight", "method", "method_label", "affinity"
  )
  group_columns <- intersect(group_columns, names(metrics))
  key_frame <- metrics[group_columns]
  key_frame[] <- lapply(key_frame, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  })
  keys <- apply(key_frame, 1L, paste, collapse = "\u001f")
  groups <- split(seq_len(nrow(metrics)), keys)

  bind_rows(lapply(groups, function(rows) {
    values <- metrics$ARI[rows]
    successful <- is.finite(values)
    output <- metrics[rows[[1L]], group_columns, drop = FALSE]
    output$attempted_replicates <- length(unique(metrics$replicate[rows]))
    output$successful_results <- sum(successful)
    output$mean_ARI <- if (any(successful)) mean(values[successful]) else NA_real_
    output$median_ARI <- if (any(successful)) stats::median(values[successful]) else NA_real_
    output$sd_ARI <- if (sum(successful) > 1L) stats::sd(values[successful]) else NA_real_
    output$mcse_mean_ARI <- if (sum(successful) > 1L) {
      output$sd_ARI / sqrt(sum(successful))
    } else {
      NA_real_
    }
    output
  }))
}

summary <- summarise_metrics(metrics)
expected_ids <- seq_along(tasks)
missing <- data.frame(task_id = setdiff(expected_ids, task_ids))
pending_methods <- bind_rows(lapply(names(design$methods), function(name) {
  specification <- design$methods[[name]]
  if (isTRUE(specification$enabled)) return(NULL)
  data.frame(
    method = name,
    label = specification$label,
    dependency_status = specification$dependency_status %||% "disabled",
    blocker = specification$blocker %||% "Disabled in manifest.",
    stringsAsFactors = FALSE
  )
}))

write.csv(full_task_table(tasks), file.path(output_dir, "design_tasks.csv"), row.names = FALSE)
write.csv(parameters, file.path(output_dir, "task_parameters.csv"), row.names = FALSE)
write.csv(status, file.path(output_dir, "task_status.csv"), row.names = FALSE)
write.csv(metrics, file.path(output_dir, "metrics.csv"), row.names = FALSE)
write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_dir, "interaction_diagnostics.csv"), row.names = FALSE)
write.csv(method_failures, file.path(output_dir, "method_failures.csv"), row.names = FALSE)
write.csv(missing, file.path(output_dir, "missing_tasks.csv"), row.names = FALSE)
write.csv(pending_methods, file.path(output_dir, "pending_methods.csv"), row.names = FALSE)

cat("Merged ", length(artifacts), " task artifact(s).\n", sep = "")
cat("Expected tasks: ", length(tasks), "\n", sep = "")
cat("Missing tasks: ", nrow(missing), "\n", sep = "")
cat("Method-level failures: ", nrow(method_failures), "\n", sep = "")
cat("Output: ", output_dir, "\n", sep = "")
