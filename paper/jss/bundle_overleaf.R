script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("bundle_overleaf.R must be run with Rscript.", call. = FALSE)
}

script_file <- sub("^--file=", "", script_argument)
paper_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
source(file.path(paper_dir, "R", "paths.R"))

ensure_output_dirs(paper_dir)
bundle_dir <- paper_path("dist", "overleaf", root = paper_dir)
rendered_candidates <- c(
  paper_path("article.tex", root = paper_dir),
  paper_path("dist", "article.tex", root = paper_dir)
)
rendered_tex <- rendered_candidates[file.exists(rendered_candidates)][1]

if (length(rendered_tex) == 0L || is.na(rendered_tex)) {
  stop(
    "Missing rendered TeX file. Checked: ",
    paste(rendered_candidates, collapse = ", "),
    ". Render article.qmd before creating the Overleaf bundle.",
    call. = FALSE
  )
}

files <- c(
  rendered_tex,
  paper_path("manydist_refs.bib", root = paper_dir),
  list.files(
    paper_path("generated", "figures", root = paper_dir),
    pattern = "[.]pdf$",
    full.names = TRUE
  )
)

copied <- file.copy(files, bundle_dir, overwrite = TRUE)
if (!all(copied)) {
  stop(
    "Could not copy every Overleaf bundle file: ",
    paste(basename(files[!copied]), collapse = ", "),
    call. = FALSE
  )
}

cat("Overleaf bundle written to ", bundle_dir, "\n", sep = "")
