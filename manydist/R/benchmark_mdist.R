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
