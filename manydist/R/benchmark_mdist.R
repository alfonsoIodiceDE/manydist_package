#' Benchmark `mdist()` over multiple method specifications
#'
#' Applies [mdist()] repeatedly over a tibble of distance-method
#' specifications, typically generated with [all_dist_method_specs()].
#'
#' Each row of `specs` is interpreted as one valid `mdist()` configuration.
#' Preset-based and custom component-based specifications are both supported.
#' Failed specifications are caught and returned in the output rather than
#' stopping the full benchmark.
#'
#' @param x A data frame or tibble of predictors, optionally including the
#'   response column.
#' @param response Optional response column inside `x`, supplied either
#'   unquoted or as a character string.
#' @param specs A tibble of method specifications. By default, this is generated
#'   with [all_dist_method_specs()]. It must contain the columns `spec_type`,
#'   `preset`, `method_cat`, `method_num`, and `commensurable`.
#'
#' @return A tibble containing the supplied specifications together with:
#' \describe{
#'   \item{result}{The corresponding output of [mdist()], or an error object if
#'   the specification failed.}
#'   \item{ok}{Logical indicator; `TRUE` if the run completed successfully,
#'   `FALSE` otherwise.}
#'   \item{error}{Error message for failed runs, `NA` otherwise.}
#' }
#'
#' @details
#' Preset specifications use the `preset` column and ignore `method_cat`,
#' `method_num`, and `commensurable`. Component specifications are evaluated as
#' `preset = "custom"` and use `method_cat`, `method_num`, and
#' `commensurable`.
#'
#' This function is intended for benchmarking, validation, and sensitivity
#' analyses across multiple distance specifications.
#'
#' @examples
#' if (requireNamespace("palmerpenguins", quietly = TRUE)) {
#'   data("penguins", package = "palmerpenguins")
#'
#'   penguins_small <- palmerpenguins::penguins |>
#'     dplyr::select(
#'       species, bill_length_mm, bill_depth_mm, flipper_length_mm,
#'       body_mass_g, island, sex
#'     ) |>
#'     tidyr::drop_na()
#'
#'   specs <- all_dist_method_specs(mode = "presets_only")
#'
#'   res <- benchmark_mdist(
#'     penguins_small,
#'     response = species,
#'     specs = specs
#'   )
#'
#'   res |>
#'     dplyr::select(spec_type, preset, ok, error)
#' }
#'
#' @export
benchmark_mdist <- function(x, response = NULL, specs = all_dist_method_specs()) {
  x <- tibble::as_tibble(x)
  specs <- tibble::as_tibble(specs)

  needed <- c(
    "spec_type", "preset",
    "method_cat", "method_num", "commensurable"
  )

  missing_cols <- setdiff(needed, names(specs))
  if (length(missing_cols) > 0) {
    stop(
      "`specs` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  response_name <- NULL
  response_quo <- rlang::enquo(response)

  if (!rlang::quo_is_null(response_quo)) {
    response_expr <- rlang::quo_get_expr(response_quo)

    if (rlang::is_string(response_expr)) {
      response_name <- response_expr
    } else {
      response_name <- rlang::as_name(response_expr)
    }

    if (!response_name %in% names(x)) {
      stop("`response` must name a column inside `x`.", call. = FALSE)
    }
  }

  results <- purrr::pmap(
    specs[needed],
    function(spec_type, preset, method_cat, method_num, commensurable) {
      tryCatch(
        {
          if (identical(spec_type, "preset")) {
            if (is.null(response_name)) {
              mdist(
                x = x,
                preset = preset
              )
            } else {
              mdist(
                x = x,
                response = response_name,
                preset = preset
              )
            }
          } else {
            if (is.null(response_name)) {
              mdist(
                x = x,
                preset = "custom",
                method_cat = method_cat,
                method_num = method_num,
                commensurable = commensurable
              )
            } else {
              mdist(
                x = x,
                response = response_name,
                preset = "custom",
                method_cat = method_cat,
                method_num = method_num,
                commensurable = commensurable
              )
            }
          }
        },
        error = function(e) e
      )
    }
  )

  specs |>
    dplyr::mutate(
      result = results,
      ok = !purrr::map_lgl(.data$result, inherits, what = "error"),
      error = purrr::map_chr(
        .data$result,
        ~ if (inherits(.x, "error")) conditionMessage(.x) else NA_character_
      )
    )
}
