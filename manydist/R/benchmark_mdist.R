#' Benchmark `mdist()` over multiple method specifications
#'
#' Apply [mdist()] repeatedly over a tibble of distance-method
#' specifications, typically generated with [all_dist_method_specs()].
#'
#' Each row of `specs` is interpreted as one valid `mdist()`
#' configuration. Preset-based and custom component-based
#' specifications are both supported.
#'
#' @param x A data frame or tibble of predictors, optionally including
#'   the response column.
#' @param response Optional response column inside `x`, supplied either
#'   unquoted or as a character string.
#' @param specs A tibble of method specifications. By default this is
#'   generated with [all_dist_method_specs()].
#'
#' @return A tibble containing the supplied specifications together with:
#' \describe{
#'   \item{result}{The corresponding output of [mdist()], or an error object
#'   if the specification failed.}
#'   \item{ok}{Logical indicator; `TRUE` if the run completed successfully,
#'   `FALSE` otherwise.}
#'   \item{error}{Error message for failed runs, `NA` otherwise.}
#' }
#'
#' @details
#' This function is intended for benchmarking, validation, and
#' sensitivity analyses across multiple distance specifications.
#'
#' @examples
#' \dontrun{
#' specs <- all_dist_method_specs()
#' res <- benchmark_mdist(Dartpoints::df, response = Name, specs = specs)
#'
#' res |>
#'   dplyr::select(spec_type, preset, distance_cont, distance_cat, ok, error)
#' }
#'
#' @export
benchmark_mdist <- function(x, response = NULL, specs = all_dist_method_specs()) {
  x <- tibble::as_tibble(x)
  specs <- tibble::as_tibble(specs)

  needed <- c(
    "spec_type", "preset", "distance_cont",
    "distance_cat", "scaling_cont", "commensurable"
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
    if (rlang::is_string(rlang::eval_tidy(response_quo))) {
      response_name <- rlang::eval_tidy(response_quo)
    } else {
      response_name <- rlang::as_name(rlang::ensym(response))
    }

    if (!response_name %in% names(x)) {
      stop("`response` must name a column inside `x`.", call. = FALSE)
    }
  }

  results <- purrr::pmap(
    specs[needed],
    function(spec_type, preset, distance_cont, distance_cat, scaling_cont, commensurable) {
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
                distance_cont = distance_cont,
                distance_cat = distance_cat,
                scaling_cont = scaling_cont,
                commensurable = commensurable
              )
            } else {
              mdist(
                x = x,
                response = response_name,
                preset = "custom",
                distance_cont = distance_cont,
                distance_cat = distance_cat,
                scaling_cont = scaling_cont,
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
