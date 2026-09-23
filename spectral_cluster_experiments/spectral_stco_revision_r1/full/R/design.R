`%||%` <- function(x, y) if (is.null(x)) y else x

full_read_design <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The full experiment runner requires the `yaml` package.", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  design <- yaml::read_yaml(path)
  if (!identical(as.integer(design$schema_version), 1L)) {
    stop("Unsupported design schema version.", call. = FALSE)
  }
  if (is.null(design$scenarios) || !length(design$scenarios)) {
    stop("The design manifest contains no scenarios.", call. = FALSE)
  }
  gate <- design$distance$categorical$interaction_gate
  if (is.null(gate) || !gate$method %in% c("permutation_max", "none")) {
    stop("The design manifest has an invalid interaction gate.", call. = FALSE)
  }
  if (!is.finite(gate$permutations) || gate$permutations < 1L ||
      gate$permutations != as.integer(gate$permutations)) {
    stop("The interaction-gate permutation count must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(gate$alpha) || gate$alpha <= 0 || gate$alpha >= 1) {
    stop("The interaction-gate alpha must lie in (0, 1).", call. = FALSE)
  }
  if (identical(gate$method, "permutation_max") &&
      1 / (as.integer(gate$permutations) + 1) > gate$alpha) {
    stop("The permutation count cannot attain the configured gate alpha.", call. = FALSE)
  }
  if (!is.finite(gate$seed_offset) || gate$seed_offset < 0 ||
      gate$seed_offset != as.integer(gate$seed_offset)) {
    stop("The interaction-gate seed offset must be a nonnegative integer.", call. = FALSE)
  }
  attr(design, "config_path") <- path
  attr(design, "config_md5") <- unname(tools::md5sum(path))
  design
}

full_expand_grid <- function(grid) {
  if (is.null(grid) || !length(grid)) return(data.frame(.row = 1L)[, FALSE])
  values <- lapply(grid, function(x) {
    if (is.null(x)) {
      stop("Grid values must be explicit; quote the string 'null' if intended.", call. = FALSE)
    }
    if (is.list(x) && !is.data.frame(x)) {
      unlist(x, recursive = FALSE, use.names = FALSE)
    } else {
      x
    }
  })
  names(values) <- names(grid)
  expanded <- do.call(
    expand.grid,
    c(values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  )
  expanded[] <- lapply(expanded, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  expanded
}

full_row_as_list <- function(x, row) {
  if (!ncol(x)) return(list())
  lapply(x[row, , drop = FALSE], function(value) value[[1L]])
}

full_expand_tasks <- function(design) {
  tasks <- list()
  task_id <- 0L
  default_replicates <- as.integer(design$execution$replicates_default)
  seed_base <- as.integer(design$execution$data_seed_base)
  scenario_stride <- as.integer(design$execution$scenario_seed_stride)
  replicate_stride <- as.integer(design$execution$replicate_seed_stride)

  enabled_names <- names(design$scenarios)[vapply(
    design$scenarios,
    function(x) isTRUE(x$enabled),
    logical(1)
  )]

  for (scenario_index in seq_along(enabled_names)) {
    scenario_name <- enabled_names[[scenario_index]]
    specification <- design$scenarios[[scenario_name]]
    outer_grid <- full_expand_grid(specification$grid)
    within_grid <- full_expand_grid(specification$within_task_grid)
    replicates <- as.integer(specification$replicates %||% default_replicates)
    if (!is.finite(replicates) || replicates < 1L) {
      stop("Scenario `", scenario_name, "` has an invalid replicate count.", call. = FALSE)
    }

    for (condition_id in seq_len(max(1L, nrow(outer_grid)))) {
      condition <- full_row_as_list(outer_grid, condition_id)
      parameters <- c(specification$fixed %||% list(), condition)
      invalid_names <- is.null(names(parameters)) || any(
        is.na(names(parameters)) | !nzchar(names(parameters)) |
          names(parameters) %in% c("FALSE", "TRUE") |
          grepl("^Var[0-9]+$", names(parameters))
      )
      if (isTRUE(invalid_names)) {
        stop(
          "Scenario `", scenario_name,
          "` has invalid parameter names; quote YAML 1.1 boolean-like keys such as 'n'.",
          call. = FALSE
        )
      }

      for (replicate_id in seq_len(replicates)) {
        task_id <- task_id + 1L
        data_seed <- seed_base + scenario_index * scenario_stride +
          replicate_id * replicate_stride
        tasks[[task_id]] <- list(
          task_id = task_id,
          scenario = scenario_name,
          scenario_index = scenario_index,
          generator = specification$generator,
          method_group = specification$method_group,
          purpose = specification$purpose,
          condition_id = condition_id,
          replicate = replicate_id,
          data_seed = as.integer(data_seed),
          parameters = parameters,
          within_task_grid = within_grid
        )
      }
    }
  }
  expected_count <- design$execution$expected_task_count
  if (!is.null(expected_count) && length(tasks) != as.integer(expected_count)) {
    stop(
      "Expanded ", length(tasks), " tasks but the manifest declares ",
      as.integer(expected_count), ".",
      call. = FALSE
    )
  }
  tasks
}

full_compact_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.logical(x)) return(tolower(as.character(x)))
  paste(as.character(x), collapse = ",")
}

full_compact_list <- function(x) {
  if (!length(x)) return("")
  paste(
    paste(names(x), vapply(x, full_compact_value, character(1)), sep = "="),
    collapse = ";"
  )
}

full_task_table <- function(tasks) {
  do.call(rbind, lapply(tasks, function(task) {
    within <- if (ncol(task$within_task_grid)) {
      paste(
        apply(task$within_task_grid, 1L, function(x) {
          full_compact_list(as.list(x))
        }),
        collapse = "|"
      )
    } else {
      ""
    }
    data.frame(
      task_id = task$task_id,
      scenario = task$scenario,
      condition_id = task$condition_id,
      replicate = task$replicate,
      data_seed = task$data_seed,
      method_group = task$method_group,
      parameters = full_compact_list(task$parameters),
      within_task_grid = within,
      stringsAsFactors = FALSE
    )
  }))
}

full_resolve_output_dir <- function(design, smoke = FALSE) {
  config_dir <- dirname(attr(design, "config_path"))
  relative <- if (isTRUE(smoke)) {
    design$outputs$smoke_directory
  } else {
    design$outputs$task_directory
  }
  config_tag <- substr(attr(design, "config_md5"), 1L, 10L)
  normalizePath(file.path(config_dir, relative, config_tag), mustWork = FALSE)
}

full_resolve_merged_dir <- function(design) {
  config_dir <- dirname(attr(design, "config_path"))
  config_tag <- substr(attr(design, "config_md5"), 1L, 10L)
  normalizePath(
    file.path(config_dir, design$outputs$merged_directory, config_tag),
    mustWork = FALSE
  )
}

full_parse_cli <- function(args, mode = c("run", "merge")) {
  mode <- match.arg(mode)
  out <- list(
    config = NULL,
    task_id = NULL,
    list_tasks = FALSE,
    count_tasks = FALSE,
    smoke = FALSE,
    include_smoke = FALSE
  )
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument %in% c("--config", "--task-id")) {
      if (index == length(args)) stop(argument, " requires a value.", call. = FALSE)
      key <- sub("^--", "", argument)
      key <- gsub("-", "_", key)
      out[[key]] <- args[[index + 1L]]
      index <- index + 2L
    } else if (argument == "--list-tasks" && mode == "run") {
      out$list_tasks <- TRUE
      index <- index + 1L
    } else if (argument == "--count-tasks" && mode == "run") {
      out$count_tasks <- TRUE
      index <- index + 1L
    } else if (argument == "--smoke" && mode == "run") {
      out$smoke <- TRUE
      index <- index + 1L
    } else if (argument == "--include-smoke" && mode == "merge") {
      out$include_smoke <- TRUE
      index <- index + 1L
    } else {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
  }
  out
}
