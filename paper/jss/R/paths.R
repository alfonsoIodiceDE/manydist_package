paper_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    candidates <- unique(c(
      current,
      file.path(current, "paper", "jss")
    ))

    for (candidate in candidates) {
      if (file.exists(file.path(candidate, "article.qmd")) &&
          file.exists(file.path(candidate, "replicate.R"))) {
        return(normalizePath(candidate, mustWork = TRUE))
      }
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  stop(
    "Could not locate paper/jss. Start inside the manydist repository ",
    "or set the working directory to paper/jss.",
    call. = FALSE
  )
}

repo_root <- function(root = paper_root()) {
  normalizePath(file.path(root, "..", ".."), mustWork = TRUE)
}

paper_path <- function(..., root = paper_root()) {
  file.path(root, ...)
}

ensure_output_dirs <- function(root = paper_root()) {
  dirs <- c(
    paper_path("generated", "figures", root = root),
    paper_path("generated", "tables", root = root),
    paper_path("generated", "results", root = root),
    paper_path("dist", "overleaf", root = root)
  )

  invisible(vapply(
    dirs,
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

