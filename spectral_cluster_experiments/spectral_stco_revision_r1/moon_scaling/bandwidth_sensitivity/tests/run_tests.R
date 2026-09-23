#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
test_dir <- dirname(script)
bandwidth_dir <- dirname(test_dir)
suite_dir <- dirname(bandwidth_dir)
revision_dir <- dirname(suite_dir)

source(file.path(revision_dir, "full", "R", "engine.R"))
source(file.path(bandwidth_dir, "R", "gaussian_affinity.R"))
testthat::test_dir(test_dir, reporter = "summary")

