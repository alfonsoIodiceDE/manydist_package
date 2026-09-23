#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
script <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "spectral_cluster_experiments/revision_r1/full/tests/run_tests.R",
    mustWork = TRUE
  )
}
full_dir <- dirname(dirname(script))
source(file.path(full_dir, "R", "design.R"))
source(file.path(full_dir, "R", "scenarios.R"))
source(file.path(full_dir, "R", "engine.R"))
testthat::test_file(file.path(full_dir, "tests", "test_design.R"), reporter = "summary")
