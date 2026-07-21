paper_distance_specifications <- function() {
  list(
    gower = list(
      preset = "gower"
    ),
    gower_commensurable = list(
      preset = "custom",
      method_cat = "matching",
      method_num = "range",
      commensurable = TRUE
    ),
    u_dep = list(
      preset = "u_dep"
    ),
    u_indep = list(
      preset = "u_indep"
    ),
    euclidean = list(
      preset = "euclidean"
    )
  )
}

distance_specifications_table <- function(specifications) {
  purrr::imap_dfr(
    specifications,
    function(specification, name) {
      tibble::tibble(
        name = name,
        preset = purrr::pluck(
          specification,
          "preset",
          .default = NA_character_
        ),
        method_cat = purrr::pluck(
          specification,
          "method_cat",
          .default = NA_character_
        ),
        method_num = purrr::pluck(
          specification,
          "method_num",
          .default = NA_character_
        ),
        commensurable = purrr::pluck(
          specification,
          "commensurable",
          .default = NA
        )
      )
    }
  )
}

validate_distance_specifications <- function(specifications) {
  mdist_arguments <- names(formals(manydist::mdist))
  registry <- manydist::dist_methods_tbl()
  registry_values <- function(argument_name) {
    registry |>
      dplyr::filter(.data$argument == .env$argument_name) |>
      dplyr::pull(.data$method)
  }

  valid_presets <- registry_values("preset")
  valid_categorical_methods <- registry_values("method_cat")
  valid_numerical_methods <- registry_values("method_num")

  purrr::iwalk(specifications, function(specification, name) {
    unknown_arguments <- setdiff(names(specification), mdist_arguments)

    if (length(unknown_arguments) > 0L) {
      stop(
        "Specification ", name, " uses unknown mdist arguments: ",
        paste(unknown_arguments, collapse = ", "),
        call. = FALSE
      )
    }

    preset <- specification$preset
    if (!is.null(preset) && !preset %in% valid_presets) {
      stop("Unknown preset in specification ", name, ": ", preset, call. = FALSE)
    }

    method_cat <- specification$method_cat
    if (!is.null(method_cat) &&
        !method_cat %in% valid_categorical_methods) {
      stop(
        "Unknown categorical method in specification ",
        name,
        ": ",
        method_cat,
        call. = FALSE
      )
    }

    method_num <- specification$method_num
    if (!is.null(method_num) &&
        !method_num %in% valid_numerical_methods) {
      stop(
        "Unknown numerical method in specification ",
        name,
        ": ",
        method_num,
        call. = FALSE
      )
    }
  })

  invisible(TRUE)
}
