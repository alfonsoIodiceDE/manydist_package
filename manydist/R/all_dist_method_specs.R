#' Create a grid of `mdist()` method specifications
#'
#' Construct a tibble of valid method specifications that can be passed to
#' [benchmark_mdist()] or iterated over manually to benchmark [mdist()] across
#' multiple configurations.
#'
#' Select named presets or construct custom specifications from the categorical
#' methods and numerical preprocessing options listed in [dist_methods_tbl()].
#'
#' @param method_cat Optional character vector restricting the categorical
#'   methods used in custom specifications.
#' @param method_num Optional character vector restricting the numerical
#'   preprocessing methods used in custom specifications.
#' @param preset `NULL` (the default) returns custom combinations only. A
#'   character vector selects named presets; `"all"` selects all named presets.
#'   Include `"custom"` explicitly to also return custom combinations, for
#'   example `c("all", "custom")`. Component filters apply only to custom rows.
#' @param commensurable Optional logical vector restricting the commensurable
#'   values used in custom specifications.
#'
#' @return A tibble where each row represents one valid `mdist()` specification.
#'   The tibble contains `preset`, `method_cat`, `method_num`, and
#'   `commensurable`.
#'
#' @details
#' Rows with `preset = "custom"` use the explicit categorical, numerical, and
#' commensurability options. Named presets determine those settings internally.
#'
#' With `commensurable = TRUE`, custom grids retain only `"none"` and
#' `"pc_scores"` numerical preprocessing. Variable-wise affine scaling cancels
#' under commensurability, so `"std"` and `"robust"` are omitted as redundant
#' with `"none"`. Range scaling is also omitted from commensurable grids; it
#' clips new observations outside the training range, so its equivalence on
#' training data does not extend to all new observations. These restrictions
#' affect grid construction, not the options accepted by [mdist()].
#'
#' Named preset rows retain `NA` in their component columns: their settings are
#' determined by the preset. Supplying component filters without requesting
#' custom rows gives an error rather than silently ignoring those filters.
#'
#' @examples
#' all_dist_method_specs()
#' all_dist_method_specs(preset = c("gower", "u_indep", "u_dep"))
#' all_dist_method_specs(preset = "all")
#' all_dist_method_specs(preset = c("all", "custom"))
#' all_dist_method_specs(method_cat = c("tvd", "le_and_ho"))
#'
#' @export
all_dist_method_specs <- function(
    method_cat = NULL,
    method_num = NULL,
    preset = NULL,
    commensurable = NULL
) {
  tbl <- dist_methods_tbl()

  preset_tbl <- tbl |>
    dplyr::filter(.data$argument == "preset", .data$method != "custom")

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

  validate_names <- function(value, choices, argument) {
    if (!is.null(value) && (!is.character(value) || anyNA(value) ||
                            any(!value %in% choices))) {
      stop("`", argument, "` must be NULL or a character vector containing: ",
           paste(choices, collapse = ", "), ".", call. = FALSE)
    }
  }
  validate_names(preset, c("all", "custom", preset_tbl$method), "preset")
  validate_names(method_cat, cat_tbl$method, "method_cat")
  validate_names(method_num, num_tbl$method, "method_num")
  if (!is.null(commensurable) &&
      (!is.logical(commensurable) || anyNA(commensurable))) {
    stop("`commensurable` must be NULL or a logical vector without missing values.",
         call. = FALSE)
  }

  include_custom <- is.null(preset) || "custom" %in% preset
  component_filters <- !is.null(method_cat) || !is.null(method_num) ||
    !is.null(commensurable)
  if (!include_custom && component_filters) {
    stop("Component filters apply to custom rows. Include 'custom' in `preset` ",
         "or leave `preset = NULL`.", call. = FALSE)
  }

  if (!"all" %in% preset) {
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
      preset = .data$method,
      method_cat = NA_character_,
      method_num = NA_character_,
      commensurable = NA
    )

  if (!include_custom) {
    return(preset_specs)
  }

  custom_specs <- tidyr::crossing(
    method_cat = cat_tbl$method,
    method_num = num_tbl$method,
    commensurable = commensurable_vals
  ) |>
    dplyr::mutate(
      preset = "custom"
    ) |>
    dplyr::select(
      preset,
      method_cat,
      method_num,
      commensurable
    ) |>
    dplyr::filter(
      !(
        .data$commensurable == TRUE &
          !.data$method_num %in% c("none", "pc_scores")
      )
    )

  dplyr::bind_rows(preset_specs, custom_specs)
}
