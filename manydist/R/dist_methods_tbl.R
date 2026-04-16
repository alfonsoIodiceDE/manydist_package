dist_methods_tbl <- function() {
  phil_methods <- philentropy::getDistMethods()

  cont_tbl <- tibble::tibble(
    method = c("euclidean", "manhattan"),
    argument = "distance_cont",
    data_type = "numeric",
    distance_basis = "numeric",
    response_aware = FALSE,
    engine = "manydist"
  )

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
    argument = "distance_cat",
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

  scaling_tbl <- tibble::tibble(
    method = c("none", "std", "pc_scores", "robust","range"),
    argument = "scaling_cont",
    data_type = "numeric",
    distance_basis = "numeric_scaling",
    response_aware = FALSE,
    engine = "manydist"
  )

  preset_tbl <- tibble::tibble(
    method = c(
      "custom", "gower", "gudmm", "dkss", "mod_gower",
      "euclidean_onehot", "hl", "u_dep", "u_indep", "u_mix"
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
      FALSE,  # euclidean_onehot
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
      "manydist",      # euclidean_onehot
      "manydist",      # hl
      "manydist",      # u_dep
      "manydist",      # u_indep
      "manydist"       # u_mix
    )
  )

  dplyr::bind_rows(cont_tbl, cat_tbl, scaling_tbl, preset_tbl) |>
    dplyr::distinct(method, argument, .keep_all = TRUE) |>
    dplyr::arrange(
      factor(argument, levels = c("distance_cont", "distance_cat", "scaling_cont", "preset")),
      factor(data_type, levels = c("numeric", "categorical", "mixed")),
      factor(
        distance_basis,
        levels = c(
          "numeric",
          "association",
          "independence",
          "indicator_scaling",
          "numeric_scaling",
          "mixed_block",
          "user_defined"
        )
      ),
      factor(engine, levels = c("manydist", "philentropy", "kdml", "external_port")),
      method
    )
}

dist_methods_tbl_cat <- function() {
  dist_methods_tbl() |>
    dplyr::filter(.data$data_type == "categorical")
}

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
