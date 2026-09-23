#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
suite_dir <- dirname(dirname(script))
revision_dir <- dirname(suite_dir)
repo_root <- dirname(dirname(revision_dir))

for (package in c("devtools", "testthat")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Tests require `", package, "`.", call. = FALSE)
  }
}
devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)
source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))
source(file.path(suite_dir, "R", "generator.R"))
source(file.path(suite_dir, "R", "gated_distance.R"))

testthat::test_file(file.path(suite_dir, "tests", "test_sim3.R"), reporter = "summary")
