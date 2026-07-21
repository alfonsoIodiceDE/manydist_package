script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("replicate.R must be run with Rscript.", call. = FALSE)
}

script_file <- sub("^--file=", "", script_argument)
paper_dir <- dirname(normalizePath(script_file, mustWork = TRUE))

source(file.path(paper_dir, "R", "paths.R"))
source(file.path(paper_dir, "R", "data.R"))
source(file.path(paper_dir, "R", "specifications.R"))
source(file.path(paper_dir, "R", "analysis.R"))
source(file.path(paper_dir, "R", "tables.R"))
source(file.path(paper_dir, "R", "replication.R"))

manifest <- run_replication(paper_dir)
cat(normalizePath(manifest, mustWork = TRUE), "\n")

