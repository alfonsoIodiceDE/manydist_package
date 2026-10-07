qmd <- normalizePath("R1_real_data.qmd", mustWork = TRUE)
lines <- readLines(qmd, warn = FALSE)
lines <- lines[seq_len(match("# Optional execution", lines) - 1L)]

chunks <- list()
inside <- FALSE
current <- character()
for (line in lines) {
  if (!inside && grepl("^```\\{r", line)) {
    inside <- TRUE
    current <- character()
  } else if (inside && grepl("^```[[:space:]]*$", line)) {
    chunks[[length(chunks) + 1L]] <- current
    inside <- FALSE
  } else if (inside && !grepl("^#\\|", line)) {
    current <- c(current, line)
  }
}
for (chunk in chunks) {
  eval(parse(text = paste(chunk, collapse = "\n")), envir = .GlobalEnv)
}

artifact_file <- r1rd_run_dataset(
  "Dermatology", match("Dermatology", dataset_registry$dataset)
)
cat(artifact_file, "\n")
