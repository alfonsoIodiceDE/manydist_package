#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "spectral_cluster_experiments/revision_r1/structured_blocks_pilot/tests/run_tests.R",
    mustWork = TRUE
  )
}
pilot_dir <- dirname(dirname(script_path))
source(file.path(pilot_dir, "R", "generator.R"))
specification <- yaml::read_yaml(file.path(pilot_dir, "design.yml"))

generated <- structured_blocks_make_pair(specification$data, 26092501L)
stopifnot(
  identical(generated$aligned$truth, generated$decoupled$truth),
  identical(
    generated$aligned$x[vapply(generated$aligned$x, is.numeric, logical(1))],
    generated$decoupled$x[vapply(generated$decoupled$x, is.numeric, logical(1))]
  ),
  nlevels(generated$aligned$truth) == 3L,
  ncol(generated$aligned$x) == 36L
)

aligned_table <- table(generated$aligned$truth, generated$aligned$categorical_regime)
decoupled_table <- table(generated$decoupled$truth, generated$decoupled$categorical_regime)
stopifnot(
  all(diag(aligned_table) == 120L),
  sum(aligned_table) == 360L,
  all(decoupled_table == 40L)
)

categorical_names <- names(generated$aligned$x)[
  vapply(generated$aligned$x, is.factor, logical(1))
]
aligned_rows <- apply(generated$aligned$x[categorical_names], 1L, paste, collapse = "|")
decoupled_rows <- apply(generated$decoupled$x[categorical_names], 1L, paste, collapse = "|")
stopifnot(identical(sort(aligned_rows), sort(decoupled_rows)))

expected_levels <- as.integer(unlist(specification$data$categorical_levels))
observed_levels <- vapply(
  generated$aligned$x[categorical_names],
  nlevels,
  integer(1)
)
stopifnot(identical(unname(observed_levels), expected_levels))

cat("Structured-block generator tests passed.\n")
