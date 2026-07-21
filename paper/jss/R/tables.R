write_csv_artifact <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path, na = "")
  normalizePath(path, mustWork = TRUE)
}

categorical_method_table <- function(root = paper_root()) {
  metadata_path <- paper_path("data", "method-metadata.csv", root = root)
  metadata <- readr::read_csv(
    metadata_path,
    show_col_types = FALSE
  ) |>
    dplyr::mutate(.order = dplyr::row_number())

  registry <- manydist::dist_methods_tbl() |>
    dplyr::filter(.data$argument == "method_cat") |>
    dplyr::select(dplyr::all_of(c(
      "method",
      "data_type",
      "distance_basis",
      "response_aware",
      "engine"
    )))

  output <- metadata |>
    dplyr::left_join(registry, by = "method") |>
    dplyr::arrange(.data$.order) |>
    dplyr::select(-dplyr::all_of(".order"))

  missing_methods <- output |>
    dplyr::filter(is.na(.data$engine)) |>
    dplyr::pull(.data$method)
  if (length(missing_methods) > 0L) {
    stop(
      "Paper method metadata contains methods absent from the package registry: ",
      paste(missing_methods, collapse = ", "),
      call. = FALSE
    )
  }

  output
}

write_replication_manifest <- function(root = paper_root()) {
  generated <- paper_path("generated", root = root)
  paths <- list.files(
    generated,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  paths <- paths[basename(paths) != "manifest.csv"]

  relative <- substring(
    normalizePath(paths, mustWork = TRUE),
    nchar(normalizePath(root, mustWork = TRUE)) + 2L
  )

  manifest <- tibble::tibble(
    path = relative,
    bytes = unname(file.info(paths)$size),
    md5 = unname(tools::md5sum(paths))
  )

  write_csv_artifact(
    manifest,
    paper_path("generated", "manifest.csv", root = root)
  )
}
