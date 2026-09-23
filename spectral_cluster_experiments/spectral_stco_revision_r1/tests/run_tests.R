args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "tests/run_tests.R"
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)

source(file.path(root, "R", "manydist_bridge.R"))
source(file.path(root, "R", "distance_r1.R"))

testthat::test_file(file.path(root, "tests", "test_distance_r1.R"), reporter = "summary")
