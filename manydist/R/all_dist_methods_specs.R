#' Create a grid of `mdist()` method specifications
#'
#' Construct a tibble of valid method specifications that can be passed to
#' [benchmark_mdist()] or iterated over manually to benchmark [mdist()] across
#' multiple configurations.
#'
#' The resulting tibble combines preset-based specifications and custom
#' component-based specifications built from the currently available
#' categorical methods and numerical preprocessing options listed in
#' [dist_methods_tbl()].
#'
#' @param mode Character string controlling the initial pool of specifications.
#'   Supported values are `"full"`, `"presets_only"`, and
#'   `"response_aware_only"`.
#' @param method_cat Optional character vector restricting the categorical
#'   methods used in component-based specifications.
#' @param method_num Optional character vector restricting the numerical
#'   preprocessing methods used in component-based specifications.
#' @param preset Optional character vector restricting preset-based
#'   specifications.
#' @param commensurable Optional logical vector restricting the commensurable
#'   values used in component-based specifications.
#'
#' @return A tibble where each row represents one valid `mdist()` specification.
#'   The tibble contains `spec_type`, `preset`, `method_cat`, `method_num`, and
#'   `commensurable`.
#'
#' @details
#' `mode` defines the initial candidate pool. Any explicit argument filters
#' supplied by the user are applied afterwards and therefore restrict the
#' selected pool further.
#'
#' With `mode = "response_aware_only"`, preset specifications and categorical
#' methods are restricted to response-aware methods.
#'
#' @examples
#' all_dist_method_specs()
#' all_dist_method_specs(mode = "presets_only")
#' all_dist_method_specs(mode = "response_aware_only")
#' all_dist_method_specs(mode = "full", method_cat = c("tvd", "le_and_ho"))
#'
#' @export
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

  num_tbl <- tbl |>
    dplyr::filter(.data$argument == "method_num")

  philentropy_to_remove <- c(
    "euclidean",
    "manhattan",
    "gower",
    "minkowski",
    "avg"
  )

  cat_tbl <- cat_tbl |>
    dplyr::filter(
      !(
        .data$engine == "philentropy" &
          .data$method %in% philentropy_to_remove
      )
    )

  if (identical(mode, "response_aware_only")) {
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

  if (identical(mode, "presets_only")) {
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
      spec_type,
      preset,
      method_cat,
      method_num,
      commensurable
    ) |>
    dplyr::filter(
      !(
        .data$commensurable == TRUE &
          !.data$method_num %in% c("std", "robust", "pc_scores")
      )
    )

  dplyr::bind_rows(preset_specs, component_specs)
}
