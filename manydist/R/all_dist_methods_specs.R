all_dist_method_specs <- function(
    mode = c("full", "presets_only", "response_aware_only"),
    method_cat = NULL,
    method_num = NULL,
    preset = NULL,
    commensurable = NULL
) {
  mode <- match.arg(mode)

  tbl <- dist_methods_tbl()

  preset_tbl <- tbl |>
    dplyr::filter(.data$argument == "preset")

  cat_tbl <- tbl |>
    dplyr::filter(.data$argument == "method_cat")

  philentropy_to_remove <- c(
    "euclidean",
    "manhattan",
    "gower",
    "minkowski",
    "avg"
  )

  cat_tbl <- cat_tbl |>
    dplyr::filter(
      !(.data$engine == "philentropy" &
          .data$method %in% philentropy_to_remove)
    )

  num_tbl <- tbl |>
    dplyr::filter(.data$argument == "method_num")

  if (mode == "response_aware_only") {
    preset_tbl <- preset_tbl |>
      dplyr::filter(.data$response_aware)

    cat_tbl <- cat_tbl |>
      dplyr::filter(.data$response_aware)
  }

  if (!is.null(preset)) {
    preset_tbl <- preset_tbl |>
      dplyr::filter(.data$method %in% preset)
  }

  if (!is.null(method_cat)) {
    cat_tbl <- cat_tbl |>
      dplyr::filter(.data$method %in% method_cat)
  }

  if (!is.null(method_num)) {
    num_tbl <- num_tbl |>
      dplyr::filter(.data$method %in% method_num)
  }

  if (is.null(commensurable)) {
    commensurable_vals <- c(TRUE, FALSE)
  } else {
    commensurable_vals <- commensurable
  }

  preset_specs <- preset_tbl |>
    dplyr::transmute(
      spec_type = "preset",
      preset = .data$method,
      method_cat = NA_character_,
      method_num = NA_character_,
      commensurable = NA
    )

  if (mode == "presets_only") {
    return(preset_specs)
  }

  component_specs <- tidyr::crossing(
    method_cat = cat_tbl$method,
    method_num = num_tbl$method,
    commensurable = commensurable_vals
  ) |>
    dplyr::mutate(
      spec_type = "component",
      preset = "custom"
    ) |>
    dplyr::select(
      .data$spec_type,
      .data$preset,
      .data$method_cat,
      .data$method_num,
      .data$commensurable
    ) |>
    dplyr::filter(
      !(.data$commensurable == TRUE &
          !.data$method_num %in% c("std", "robust", "pc_scores"))
    )

  dplyr::bind_rows(preset_specs, component_specs)
}
