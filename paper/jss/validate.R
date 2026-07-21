script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("validate.R must be run with Rscript.", call. = FALSE)
}

script_file <- sub("^--file=", "", script_argument)
paper_dir <- dirname(normalizePath(script_file, mustWork = TRUE))

r_source_files <- sort(list.files(
  file.path(paper_dir, "R"),
  pattern = "[.]R$",
  full.names = TRUE
))

entrypoint_files <- file.path(
  paper_dir,
  c("_targets.R", "bundle_overleaf.R", "replicate.R")
)
parse_files <- c(r_source_files, entrypoint_files)

for (file in parse_files) {
  parse(file = file)
}

source(file.path(paper_dir, "R", "paths.R"))
source(file.path(paper_dir, "R", "data.R"))
source(file.path(paper_dir, "R", "specifications.R"))
source(file.path(paper_dir, "R", "tables.R"))

required_files <- c(
  "README.md",
  "MIGRATION.md",
  ".gitignore",
  "Makefile",
  "_quarto.yml",
  "_targets.R",
  "article.qmd",
  "manydist_refs.bib",
  "replicate.R",
  "bundle_overleaf.R",
  "data/related-software.csv",
  "data/method-metadata.csv"
)

missing_files <- required_files[
  !file.exists(file.path(paper_dir, required_files))
]
if (length(missing_files) > 0L) {
  stop(
    "Missing required paper files: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

load_manydist_source(paper_dir)

penguins <- prepare_penguins()
if (!identical(dim(penguins$predictors), c(333L, 6L))) {
  stop(
    "Unexpected penguin predictor dimensions: ",
    paste(dim(penguins$predictors), collapse = " x "),
    call. = FALSE
  )
}

specifications <- paper_distance_specifications()
validate_distance_specifications(specifications)
invisible(categorical_method_table(paper_dir))

cat("Validated ", length(parse_files), " R source files.\n", sep = "")
cat("Penguin predictors: 333 x 6.\n")
cat("Distance specifications: ", length(specifications), ".\n", sep = "")
cat("Paper scaffold validation passed.\n")
