#' List available `manydist` methods
#'
#' Returns a tibble describing the categorical methods, numerical preprocessing
#' options, and preset specifications currently available in `manydist`.
#'
#' The table is used by helpers such as [all_dist_method_specs()] and
#' [benchmark_mdist()] to construct valid `mdist()` specifications.
#'
#' @return A tibble with one row per available method. The columns are:
#' \describe{
#'   \item{method}{Name of the method or preset.}
#'   \item{argument}{The `mdist()` argument to which the method belongs:
#'   `"method_cat"`, `"method_num"`, or `"preset"`.}
#'   \item{data_type}{Type of variables targeted by the method.}
#'   \item{distance_basis}{Broad methodological family.}
#'   \item{response_aware}{Logical; whether the method can use a response
#'   variable.}
#'   \item{engine}{Implementation source.}
#' }
#'
#' @examples
#' dist_methods_tbl()
#'
#' dist_methods_tbl() |>
#'   dplyr::count(argument)
#'
#' @export
dist_methods_tbl <- function() {
  phil_methods <- philentropy::getDistMethods()

  cat_tbl <- tibble::tibble(
    method = c(
      # association-based custom
      "tvd", "le_and_ho", "gifi_chi2",

      # association-based via philentropy
      phil_methods,

      # independence-based
      "matching", "eskin", "goodall_3", "goodall_4",
      "iof", "of", "lin", "var_entropy", "var_mutability",

      # indicator-scaling based
      "none", "st_dev", "HL", "cat_dis", "HLeucl", "mca"
    ),
    argument = "method_cat",
    data_type = "categorical",
    distance_basis = c(
      rep("association", 3 + length(phil_methods)),
      rep("independence", 9),
      rep("indicator_scaling", 6)
    ),
    response_aware = c(
      rep(TRUE, 3 + length(phil_methods)),
      rep(FALSE, 9),
      rep(FALSE, 6)
    ),
    engine = c(
      rep("manydist", 3),
      rep("philentropy", length(phil_methods)),
      rep("manydist", 9),
      rep("manydist", 6)
    )
  )

  num_tbl <- tibble::tibble(
    method = c("none", "std", "pc_scores", "robust", "range"),
    argument = "method_num",
    data_type = "numeric",
    distance_basis = "numeric_preprocessing",
    response_aware = FALSE,
    engine = "manydist"
  )

  preset_tbl <- tibble::tibble(
    method = c(
      "custom", "gower", "gudmm", "dkss", "mod_gower",
      "euclidean", "hl", "u_dep", "u_indep", "u_mix"
    ),
    argument = "preset",
    data_type = "mixed",
    distance_basis = c(
      "user_defined",
      rep("mixed_block", 9)
    ),
    response_aware = c(
      FALSE,  # custom
      FALSE,  # gower
      FALSE,  # gudmm
      FALSE,  # dkss
      FALSE,  # mod_gower
      FALSE,  # euclidean
      FALSE,  # hl
      TRUE,   # u_dep
      FALSE,  # u_indep
      TRUE    # u_mix
    ),
    engine = c(
      "manydist",      # custom
      "manydist",      # gower
      "external_port", # gudmm
      "kdml",          # dkss
      "external_port", # mod_gower
      "manydist",      # euclidean
      "manydist",      # hl
      "manydist",      # u_dep
      "manydist",      # u_indep
      "manydist"       # u_mix
    )
  )

  dplyr::bind_rows(cat_tbl, num_tbl, preset_tbl) |>
    dplyr::distinct(.data$method, .data$argument, .keep_all = TRUE) |>
    dplyr::arrange(
      factor(.data$argument, levels = c("method_cat", "method_num", "preset")),
      factor(.data$data_type, levels = c("numeric", "categorical", "mixed")),
      factor(
        .data$distance_basis,
        levels = c(
          "association",
          "independence",
          "indicator_scaling",
          "numeric_preprocessing",
          "mixed_block",
          "user_defined"
        )
      ),
      factor(.data$engine, levels = c("manydist", "philentropy", "kdml", "external_port")),
      .data$method
    )
}


#' List available categorical dissimilarities
#'
#' Convenience wrapper around [dist_methods_tbl()] returning only methods that
#' can be used as categorical dissimilarities through `method_cat`.
#'
#' @return A tibble containing the rows of [dist_methods_tbl()] with
#'   `argument == "method_cat"`.
#'
#' @examples
#' dist_methods_tbl_cat()
#'
#' @export
dist_methods_tbl_cat <- function() {
  dist_methods_tbl() |>
    dplyr::filter(.data$argument == "method_cat")
}


#' List response-aware methods
#'
#' Returns the methods in [dist_methods_tbl()] that can use a response variable.
#' The output can optionally be restricted by data type or by `mdist()`
#' argument.
#'
#' @param data_type Optional character vector used to restrict the returned
#'   methods by data type, for example `"categorical"` or `"mixed"`.
#' @param argument Optional character vector used to restrict the returned
#'   methods by argument, for example `"method_cat"` or `"preset"`.
#'
#' @return A character vector of response-aware method names.
#'
#' @examples
#' response_aware_methods()
#' response_aware_methods(argument = "method_cat")
#' response_aware_methods(argument = "preset")
#'
#' @export
response_aware_methods <- function(data_type = NULL, argument = NULL) {
  out <- dist_methods_tbl()

  if (!is.null(data_type)) {
    out <- out |>
      dplyr::filter(.data$data_type %in% data_type)
  }

  if (!is.null(argument)) {
    out <- out |>
      dplyr::filter(.data$argument %in% argument)
  }

  out |>
    dplyr::filter(.data$response_aware) |>
    dplyr::pull(.data$method)
}
