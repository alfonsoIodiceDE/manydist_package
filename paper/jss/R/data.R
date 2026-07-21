load_manydist_source <- function(root = paper_root()) {
  package_dir <- file.path(repo_root(root), "manydist")

  if (requireNamespace("pkgload", quietly = TRUE) &&
      file.exists(file.path(package_dir, "DESCRIPTION"))) {
    pkgload::load_all(package_dir, quiet = TRUE)
  } else if (!requireNamespace("manydist", quietly = TRUE)) {
    stop(
      "The manydist package is not installed and pkgload is unavailable.",
      call. = FALSE
    )
  }

  invisible(as.character(utils::packageVersion("manydist")))
}

prepare_penguins <- function() {
  if (!requireNamespace("palmerpenguins", quietly = TRUE)) {
    stop("Package palmerpenguins is required for paper replication.", call. = FALSE)
  }

  columns <- c(
    "species",
    "bill_length_mm",
    "bill_depth_mm",
    "flipper_length_mm",
    "body_mass_g",
    "island",
    "sex"
  )

  data <- palmerpenguins::penguins |>
    dplyr::select(dplyr::all_of(columns)) |>
    tidyr::drop_na()

  list(
    full = data,
    predictors = data |>
      dplyr::select(-dplyr::all_of("species")),
    response = data |>
      dplyr::pull(.data$species)
  )
}

penguins_summary_table <- function(penguins) {
  purrr::imap_dfr(
    penguins$full,
    function(value, variable) {
      numeric_variable <- is.numeric(value)

      tibble::tibble(
        variable = variable,
        type = if (numeric_variable) "numeric" else "categorical",
        n = length(value),
        levels = if (is.factor(value)) nlevels(value) else NA_integer_,
        mean = if (numeric_variable) mean(value) else NA_real_,
        sd = if (numeric_variable) stats::sd(value) else NA_real_
      )
    }
  )
}

paper_seeds <- function() {
  c(
    distance = 20260721L,
    lovo = 20260722L,
    knn_split = 20260723L,
    knn_tuning = 20260724L,
    spectral = 20260725L
  )
}
