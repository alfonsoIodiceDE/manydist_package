#!/usr/bin/env Rscript

full_script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(argument)) {
    return(normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE))
  }
  normalizePath(
    "spectral_cluster_experiments/revision_r1/full/run_task.R",
    mustWork = TRUE
  )
}

script_file <- full_script_path()
full_dir <- dirname(script_file)
revision_dir <- dirname(full_dir)
repo_root <- dirname(dirname(revision_dir))

source(file.path(full_dir, "R", "design.R"))

arguments <- full_parse_cli(commandArgs(trailingOnly = TRUE), mode = "run")
config_path <- arguments$config %||% file.path(full_dir, "design.yml")
design <- full_read_design(config_path)
tasks <- full_expand_tasks(design)

if (isTRUE(arguments$count_tasks)) {
  cat(length(tasks), "\n", sep = "")
  quit(status = 0L)
}

if (isTRUE(arguments$list_tasks)) {
  write.table(
    full_task_table(tasks),
    file = stdout(),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  quit(status = 0L)
}

if (is.null(arguments$task_id)) {
  stop("Supply exactly one `--task-id`, or use `--list-tasks`/`--count-tasks`.", call. = FALSE)
}
task_id <- suppressWarnings(as.integer(arguments$task_id))
if (!is.finite(task_id) || task_id < 1L || task_id > length(tasks)) {
  stop("`--task-id` is outside the manifest task range.", call. = FALSE)
}
task <- tasks[[task_id]]

output_dir <- full_resolve_output_dir(design, smoke = arguments$smoke)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file_prefix <- if (isTRUE(arguments$smoke)) "smoke-task" else "task"
output_file <- file.path(output_dir, sprintf("%s-%05d.rds", file_prefix, task_id))

if (file.exists(output_file)) {
  cat("Task output already exists; leaving it unchanged:\n", output_file, "\n", sep = "")
  quit(status = 0L)
}

source(file.path(full_dir, "R", "scenarios.R"))

required_packages <- c("devtools", "mclust", "yaml")
enabled_scenarios <- design$scenarios[vapply(
  design$scenarios,
  function(x) isTRUE(x$enabled),
  logical(1)
)]
enabled_generators <- vapply(enabled_scenarios, `[[`, character(1), "generator")
if ("legacy_simulation_1" %in% enabled_generators) {
  required_packages <- c(required_packages, "mvtnorm")
}
enabled_engines <- vapply(
  design$methods[vapply(design$methods, function(x) isTRUE(x$enabled), logical(1))],
  `[[`,
  character(1),
  "engine"
)
if ("clustMixType_kprototypes" %in% enabled_engines) {
  required_packages <- c(required_packages, "clustMixType")
}
if ("manydist_dkss_distance" %in% enabled_engines) {
  required_packages <- c(required_packages, "kdml")
}
required_packages <- unique(required_packages)
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)
source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))
source(file.path(revision_dir, "moon_scaling", "R", "gated_distance.R"))
source(file.path(full_dir, "R", "engine.R"))

warnings_seen <- character()
started <- Sys.time()
run_error <- NA_character_
generated <- NULL
evaluated <- NULL

withCallingHandlers(
  tryCatch(
    {
      generated <- full_generate_scenario(
        task,
        design,
        smoke = isTRUE(arguments$smoke)
      )
      evaluated <- full_evaluate_task(task, generated, design, repo_root)
    },
    error = function(e) {
      run_error <<- conditionMessage(e)
    }
  ),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

finished <- Sys.time()
method_errors <- if (!is.null(evaluated) && nrow(evaluated$results)) {
  unique(stats::na.omit(evaluated$results$error))
} else {
  character()
}
status <- if (!is.na(run_error)) {
  "failed"
} else if (length(method_errors)) {
  "complete_with_method_failures"
} else {
  "complete"
}

package_versions <- vapply(
  c("manydist", "mclust", "yaml", "devtools", "clustMixType", "kdml"),
  function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  },
  character(1)
)

artifact <- list(
  schema_version = 1L,
  analysis_id = design$analysis_id,
  config_path = attr(design, "config_path"),
  config_md5 = attr(design, "config_md5"),
  smoke = isTRUE(arguments$smoke),
  status = status,
  task = task,
  effective_parameters = if (is.null(generated)) NULL else generated$parameters,
  generator_diagnostics = if (is.null(generated)) NULL else generated$diagnostics,
  results = if (is.null(evaluated)) data.frame() else evaluated$results,
  interaction_diagnostics = if (is.null(evaluated)) data.frame() else evaluated$diagnostics,
  r1_provenance = if (is.null(evaluated)) NULL else evaluated$r1_provenance,
  warnings = unique(warnings_seen),
  error = run_error,
  started_at = format(started, tz = "UTC", usetz = TRUE),
  finished_at = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
  runtime = list(
    R_version = R.version.string,
    platform = R.version$platform,
    package_versions = package_versions,
    session_info = capture.output(utils::sessionInfo())
  )
)

temporary_file <- tempfile(
  pattern = paste0(".", file_prefix, "-", sprintf("%05d", task_id), "-"),
  tmpdir = output_dir,
  fileext = ".rds"
)
saveRDS(artifact, temporary_file, version = 3)
if (file.exists(output_file)) {
  unlink(temporary_file)
  stop("Task output appeared during the run; refusing to overwrite it.", call. = FALSE)
}
if (!file.rename(temporary_file, output_file)) {
  unlink(temporary_file)
  stop("Could not atomically move the task artifact into place.", call. = FALSE)
}

cat("Task ", task_id, " ", status, "\n", sep = "")
cat("Scenario: ", task$scenario, "\n", sep = "")
cat("Output: ", output_file, "\n", sep = "")
cat("Elapsed seconds: ", sprintf("%.3f", artifact$elapsed_seconds), "\n", sep = "")

if (identical(status, "failed")) quit(status = 1L)
