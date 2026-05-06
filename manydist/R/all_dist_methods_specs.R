#' Create a grid of `mdist()` method specifications
#'
#' Construct a tibble of valid method specifications that can be passed
#' to [benchmark_mdist()] or iterated over manually to benchmark
#' [mdist()] across multiple configurations.
#'
#' The resulting tibble combines preset-based specifications and
#' component-based specifications built from the currently available
#' continuous, categorical, and scaling options listed in
#' [dist_methods_tbl()].
#'
#' @param mode Character string controlling the initial pool of
#'   specifications. Supported values are:
#'   \describe{
#'     \item{"full"}{Return all presets and all component-based combinations.}
#'     \item{"presets_only"}{Return preset-based specifications only.}
#'     \item{"response_aware_only"}{Return response-aware presets and
#'     component-based specifications using only response-aware
#'     categorical distances.}
#'   }
#' @param distance_cont Optional character vector restricting the
#'   `distance_cont` values used in component-based specifications.
#' @param distance_cat Optional character vector restricting the
#'   `distance_cat` values used in component-based specifications.
#' @param scaling_cont Optional character vector restricting the
#'   `scaling_cont` values used in component-based specifications.
#' @param preset Optional character vector restricting preset-based
#'   specifications.
#' @param commensurable Optional logical vector restricting the
#'   `commensurable` values used in component-based specifications.
#'
#' @return A tibble where each row represents one valid `mdist()`
#' specification. The tibble contains:
#' \describe{
#'   \item{spec_type}{Either `"preset"` or `"component"`.}
#'   \item{preset}{Preset name, if applicable; otherwise `NA`.}
#'   \item{distance_cont}{Continuous-distance specification, if applicable.}
#'   \item{distance_cat}{Categorical-distance specification, if applicable.}
#'   \item{scaling_cont}{Continuous scaling specification, if applicable.}
#'   \item{commensurable}{Logical flag for component-based specifications;
#'   `NA` for preset-based specifications.}
#' }
#'
#' @details
#' `mode` defines the initial candidate pool. Any explicit argument
#' filters supplied by the user are applied afterwards and therefore
#' restrict the selected pool further.
#'
#' @examples
#' all_dist_method_specs()
#' all_dist_method_specs(mode = "presets_only")
#' all_dist_method_specs(mode = "response_aware_only")
#' all_dist_method_specs(mode = "full", distance_cat = c("tvd", "le_and_ho"))
#'
#' @export
all_dist_method_specs <- function(
    mode = c("full", "presets_only", "response_aware_only"),
    distance_cont = NULL,
    distance_cat = NULL,
    scaling_cont = NULL,
    preset = NULL,
    commensurable = NULL
) {
  mode <- match.arg(mode)

  tbl <- dist_methods_tbl()

  preset_tbl <- tbl |>
    dplyr::filter(.data$argument == "preset")

  cont_tbl <- tbl |>
    dplyr::filter(.data$argument == "distance_cont")

  philentropy_to_remove <- c(
    "euclidean",
    "manhattan",
    "gower",
    "minkowski",
    "avg"
  )

  cat_tbl <- tbl |>
    dplyr::filter(.data$argument == "distance_cat",
                  !(engine == "philentropy" & method %in% philentropy_to_remove))



  scaling_tbl <- tbl |>
    dplyr::filter(.data$argument == "scaling_cont")

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

  if (!is.null(distance_cont)) {
    cont_tbl <- cont_tbl |>
      dplyr::filter(.data$method %in% distance_cont)
  }

  if (!is.null(distance_cat)) {
    cat_tbl <- cat_tbl |>
      dplyr::filter(.data$method %in% distance_cat)
  }

  if (!is.null(scaling_cont)) {
    scaling_tbl <- scaling_tbl |>
      dplyr::filter(.data$method %in% scaling_cont)
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
      distance_cont = NA_character_,
      distance_cat = NA_character_,
      scaling_cont = NA_character_,
      commensurable = NA
    )

  if (mode == "presets_only") {
    return(preset_specs)
  }

  component_specs <- tidyr::crossing(
    distance_cont = cont_tbl$method,
    distance_cat = cat_tbl$method,
    scaling_cont = scaling_tbl$method,
    commensurable = commensurable_vals
  ) |>
    dplyr::mutate(
      spec_type = "component",
      preset = "custom"
    ) |>
    dplyr::select(
      spec_type,
      preset,
      distance_cont,
      distance_cat,
      scaling_cont,
      commensurable
    ) |> dplyr::filter(
    !(commensurable == TRUE & !scaling_cont %in% c("std", "robust", "pc_scores"))
  )




  dplyr::bind_rows(preset_specs, component_specs)
}
