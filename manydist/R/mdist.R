`%||%` <- function(a, b) if (!is.null(a)) a else b


.safe_comm <- function(b_v_d) {
  den <- mean(b_v_d)
  if (!is.finite(den) || den <= .Machine$double.eps) den <- 1
  b_v_d / den
}


.to_dissimilarity <- function(dist_matrix, reference = NULL) {
  if (nrow(dist_matrix) == ncol(dist_matrix) &&
      isTRUE(all.equal(dist_matrix, t(dist_matrix), tolerance = 1e-10))) {

    d <- as.dist(dist_matrix)
    class(d) <- c("dissimilarity", "dist")

    if (!is.null(reference)) {
      for (att in intersect(names(attributes(reference)),
                            c("Labels", "Size", "Diag", "Upper"))) {
        attr(d, att) <- attr(reference, att)
      }
    }

  } else {
    d <- dist_matrix
    class(d) <- c("dissimilarity", "matrix")
  }

  d
}


.resolve_preset <- function(
    preset,
    distance_cont,
    distance_cat,
    commensurable,
    scaling_cont,
    interaction
) {

  preset_defs <- list(
    custom = list(),

    unbiased_dependent = list(
      distance_cont = "manhattan",
      distance_cat  = "tvd",
      commensurable = TRUE,
      scaling_cont  = "pc_scores",
      interaction   = FALSE
    ),
    hl = list(
      distance_cont = "euclidean",
      distance_cat  = "HLeucl",
      commensurable = FALSE,
      scaling_cont  = "std",
      ncomp         = NULL,
      threshold     = NULL
    ),
    u_dep = list(
      distance_cont = "manhattan",
      distance_cat  = "tvd",
      commensurable = TRUE,
      scaling_cont  = "pc_scores",
      interaction   = FALSE
    ),
    u_indep = list(
      distance_cont = "manhattan",
      distance_cat  = "matching",
      commensurable = TRUE,
      scaling_cont  = "std",
      interaction   = FALSE
    ),
    u_mix = list(
      distance_cont = "manhattan",
      distance_cat  = "tvd",
      commensurable = TRUE,
      scaling_cont  = "std",
      interaction   = FALSE
    )
  )

  if (!preset %in% names(preset_defs)) {
    return(list(
      distance_cont = distance_cont,
      distance_cat  = distance_cat,
      commensurable = commensurable,
      scaling_cont  = scaling_cont,
      interaction   = interaction
    ))
  }

  defs <- preset_defs[[preset]]

  list(
    distance_cont = defs$distance_cont %||% distance_cont,
    distance_cat  = defs$distance_cat  %||% distance_cat,
    commensurable = defs$commensurable %||% commensurable,
    scaling_cont  = defs$scaling_cont  %||% scaling_cont,
    interaction   = defs$interaction   %||% interaction
  )
}


.mdist_generic <- function(
    cont_data, cat_data,
    cont_data_val = NULL, cat_data_val = NULL,
    y = NULL,
    distance_cont, distance_cat,
    commensurable, scaling_cont,
    ncomp, threshold,
    interaction, prop_nn, score, decision
) {

  cont_dist_mat <- NULL
  cat_dist_mat  <- NULL

  if (!is.null(cont_data)) {
    cont_dist_mat <- ndist(
      x = cont_data,
      validate_x = cont_data_val,
      method = distance_cont,
      commensurable = commensurable,
      scaling = scaling_cont,
      ncomp = ncomp,
      threshold = threshold
    ) |>
      as.matrix()
  }

  if (!is.null(cat_data)) {
    cat_dist_mat <- cdist(
      x = cat_data,
      validate_x = cat_data_val,
      response = y,
      method = distance_cat,
      commensurable = commensurable
    )$distance_mat |>
      as.matrix()
  }

  if (!is.null(cont_dist_mat) && !is.null(cat_dist_mat)) {
    distance_mat <- cont_dist_mat + cat_dist_mat
  } else if (!is.null(cont_dist_mat)) {
    distance_mat <- cont_dist_mat
  } else if (!is.null(cat_dist_mat)) {
    distance_mat <- cat_dist_mat
  } else {
    stop("No variables available to compute distance.", call. = FALSE)
  }

  if (isTRUE(interaction)) {
    if (is.null(cont_dist_mat) || is.null(cat_dist_mat) || is.null(cat_data) || is.null(cont_data)) {
      warning("`interaction = TRUE` requires mixed data; ignored.", call. = FALSE)
    } else if (!is.null(cont_data_val) || !is.null(cat_data_val)) {
      warning("`interaction = TRUE` currently only implemented for train-train distances; ignored.", call. = FALSE)
    } else {
      i_distance_matrix <- idist(
        D = cont_dist_mat,
        cat_data = cat_data,
        pi_nn = prop_nn,
        score = score,
        decision = decision
      )

      distance_mat <- cont_dist_mat +
        (ncol(cat_data) / ncol(cont_data)) * (cat_dist_mat + i_distance_matrix)
    }
  }

  if (!is.null(cont_dist_mat) && !is.null(cat_dist_mat) &&
      identical(distance_cont, "euclidean") && identical(distance_cat, "HLeucl")) {
    distance_mat <- sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
  }

  distance_mat
}


.prep_mdist <- function(x, response = NULL,
                        distance_cont = "manhattan", distance_cat = "tvd",
                        commensurable = TRUE, scaling_cont = "std",
                        ncomp = NULL, threshold = NULL,
                        preset = "custom", interaction = FALSE,
                        prop_nn = 0.1, score = "ba", decision = "prior_corrected") {

  if (identical(distance_cat, "tot_var_dist")) {
    distance_cat <- "tvd"
  }

  resp_name <- NULL
  response_quo <- rlang::enquo(response)

  if (!rlang::quo_is_null(response_quo)) {
    resp_idx <- tidyselect::eval_select(response_quo, x)

    if (length(resp_idx) != 1L) {
      stop("`response` must select exactly one column.", call. = FALSE)
    }

    resp_name <- names(resp_idx)
    y <- x[[resp_name]]
    x <- x[, setdiff(colnames(x), resp_name), drop = FALSE]
  } else {
    y <- NULL
  }

  x <- tibble::as_tibble(x)

  if (is.null(ncomp)) ncomp <- ncol(x)

  args <- .resolve_preset(
    preset = preset,
    distance_cont = distance_cont,
    distance_cat  = distance_cat,
    commensurable = commensurable,
    scaling_cont  = scaling_cont,
    interaction   = interaction
  )

  distance_cont <- args$distance_cont
  distance_cat  <- args$distance_cat
  commensurable <- args$commensurable
  scaling_cont  <- args$scaling_cont
  interaction   <- args$interaction

  cat_data  <- x |> dplyr::select(where(is.factor))
  cont_data <- x |> dplyr::select(where(is.numeric))

  if (ncol(cat_data) == 0) cat_data <- NULL
  if (ncol(cont_data) == 0) cont_data <- NULL

  if (!is.null(cat_data) && identical(distance_cat, "tvd") && ncol(cat_data) == 1) {
    warning("'tvd' requires >1 categorical variable. Switching to 'matching'.", call. = FALSE)
    distance_cat <- "matching"
  }

  if (!is.null(cont_data) && identical(scaling_cont, "pc_scores") && ncol(cont_data) == 1) {
    warning(
      "With 1 continuous variable, 'pc_scores' is equivalent to standardization. Switching to scaling_cont='std'.",
      call. = FALSE
    )
    scaling_cont <- "std"
  }

  cat_levels <- NULL
  if (!is.null(cat_data)) {
    cat_levels <- purrr::map(cat_data, levels)
  }

  # objects fitted on training and reused later
  gower_prep <- NULL
  dummy_recipe <- NULL

  if (preset == "gower") {
    if (is.null(cont_data) && !is.null(cat_data)) {
      gower_prep <- gower_recipe(data = cat_data) |> recipes::prep(training = cat_data)
    } else if (!is.null(cont_data) && is.null(cat_data)) {
      gower_prep <- gower_recipe(data = cont_data) |> recipes::prep(training = cont_data)
    } else if (!is.null(cont_data) && !is.null(cat_data)) {
      gower_prep <- gower_recipe(data = x) |> recipes::prep(training = x)
    }
  }

  if (preset == "euclidean_onehot" && !is.null(cat_data)) {
    dummy_recipe <- recipes::recipe(~ ., data = cat_data) |>
      recipes::step_dummy(recipes::all_nominal(), one_hot = TRUE) |>
      recipes::prep(training = cat_data)
  }

  structure(
    list(
      x_train        = x,
      y              = y,
      response_col   = resp_name,
      preset         = preset,
      distance_cont  = distance_cont,
      distance_cat   = distance_cat,
      commensurable  = commensurable,
      scaling_cont   = scaling_cont,
      ncomp          = ncomp,
      threshold      = threshold,
      interaction    = interaction,
      prop_nn        = prop_nn,
      score          = score,
      decision       = decision,
      cont_data      = cont_data,
      cat_data       = cat_data,
      cont_names     = names(cont_data),
      cat_names      = names(cat_data),
      cat_levels     = cat_levels,
      gower_prep     = gower_prep,
      dummy_recipe   = dummy_recipe
    ),
    class = "mdist_preprocessor"
  )
}


.apply_mdist <- function(prep, new_data = NULL) {
  stopifnot(inherits(prep, "mdist_preprocessor"))

  validate_x <- NULL
  cont_data_val <- NULL
  cat_data_val <- NULL

  if (!is.null(new_data)) {
    validate_x <- tibble::as_tibble(new_data)

    if (!is.null(prep$response_col) && prep$response_col %in% names(validate_x)) {
      validate_x <- validate_x[, setdiff(colnames(validate_x), prep$response_col), drop = FALSE]
    }

    # keep only training columns and in training order
    validate_x <- validate_x |>
      dplyr::select(dplyr::all_of(colnames(prep$x_train)))

    # harmonize factor levels
    if (!is.null(prep$cat_names)) {
      validate_x <- validate_x |>
        dplyr::mutate(dplyr::across(
          dplyr::all_of(prep$cat_names),
          ~ factor(as.character(.x), levels = prep$cat_levels[[dplyr::cur_column()]])
        ))
    }

    if (!is.null(prep$cat_names)) {
      cat_data_val <- validate_x |>
        dplyr::select(dplyr::all_of(prep$cat_names))
      if (ncol(cat_data_val) == 0) cat_data_val <- NULL
    }

    if (!is.null(prep$cont_names)) {
      cont_data_val <- validate_x |>
        dplyr::select(dplyr::all_of(prep$cont_names))
      if (ncol(cont_data_val) == 0) cont_data_val <- NULL
    }
  }

  cont_data <- prep$cont_data
  cat_data  <- prep$cat_data

  generic_presets <- c("custom", "unbiased_dependent", "u_dep", "u_indep", "u_mix")

  if (prep$preset %in% generic_presets) {

    distance_mat <- .mdist_generic(
      cont_data = cont_data,
      cat_data = cat_data,
      cont_data_val = cont_data_val,
      cat_data_val = cat_data_val,
      y = prep$y,
      distance_cont = prep$distance_cont,
      distance_cat = prep$distance_cat,
      commensurable = prep$commensurable,
      scaling_cont = prep$scaling_cont,
      ncomp = prep$ncomp,
      threshold = prep$threshold,
      interaction = prep$interaction,
      prop_nn = prep$prop_nn,
      score = prep$score,
      decision = prep$decision
    )

  } else if (prep$preset == "gower") {

    if (is.null(validate_x)) {

      if (is.null(cont_data) && !is.null(cat_data)) {

        if (isFALSE(prep$commensurable)) {
          distance_mat <- ncol(cat_data) * cluster::daisy(cat_data, metric = "gower") |> as.matrix()
        } else {
          gowerlist <- cat_data |>
            purrr::map(~ cluster::daisy(tibble::as_tibble(.x), metric = "gower") |> as.matrix())
          gowerlist <- tibble::tibble(gowdist = gowerlist) |>
            dplyr::mutate(commgow = purrr::map(.data$gowdist, ~ .safe_comm(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)
        }

      } else if (!is.null(cont_data) && is.null(cat_data)) {

        if (isFALSE(prep$commensurable)) {
          distance_mat <- ncol(cont_data) * cluster::daisy(cont_data, metric = "gower") |> as.matrix()
        } else {
          gowerlist <- cont_data |>
            purrr::map(~ cluster::daisy(tibble::as_tibble(.x), metric = "gower") |> as.matrix())
          gowerlist <- tibble::tibble(gowdist = gowerlist) |>
            dplyr::mutate(commgow = purrr::map(.data$gowdist, ~ .safe_comm(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)
        }

      } else if (!is.null(cont_data) && !is.null(cat_data)) {

        if (isFALSE(prep$commensurable)) {
          distance_mat <- cluster::daisy(prep$x_train, metric = "gower") |> as.matrix()
        } else {
          gowerlist <- prep$x_train |>
            purrr::map(~ cluster::daisy(tibble::as_tibble(.x), metric = "gower") |> as.matrix())
          gowerlist <- tibble::tibble(gowdist = gowerlist) |>
            dplyr::mutate(commgow = purrr::map(.data$gowdist, ~ .safe_comm(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)
        }

      } else {
        stop("No variables available to compute distance.", call. = FALSE)
      }

    } else {

      train_gow <- recipes::bake(prep$gower_prep, new_data = NULL)
      valid_gow <- recipes::bake(prep$gower_prep, new_data = validate_x)

      gowerlist <- purrr::map2(
        .x = tibble::as_tibble(train_gow),
        .y = tibble::as_tibble(valid_gow),
        .f = function(x, xnew) {
          b_v_d <- abs(outer(as.numeric(xnew), as.numeric(x), "-"))
          if (isTRUE(prep$commensurable)) .safe_comm(b_v_d) else b_v_d
        }
      )
      distance_mat <- Reduce(`+`, gowerlist)
    }

  } else if (prep$preset == "euclidean_onehot") {

    distance_cont <- "euclidean"
    commensurable <- FALSE
    scaling_cont  <- "std"

    if (is.null(cat_data) && !is.null(cont_data)) {

      distance_mat <- ndist(
        x = cont_data,
        validate_x = cont_data_val,
        method = distance_cont,
        commensurable = commensurable,
        scaling = scaling_cont
      ) |>
        as.matrix()

    } else if (!is.null(cat_data) && is.null(cont_data)) {

      cat_data_dummy <- recipes::bake(prep$dummy_recipe, new_data = NULL)

      if (is.null(validate_x)) {
        distance_mat <- ndist(
          x = cat_data_dummy,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling_cont
        ) |>
          as.matrix()
      } else {
        cat_data_val_dummy <- recipes::bake(prep$dummy_recipe, new_data = cat_data_val)

        distance_mat <- ndist(
          x = cat_data_dummy,
          validate_x = cat_data_val_dummy,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling_cont
        ) |>
          as.matrix()
      }

    } else if (!is.null(cat_data) && !is.null(cont_data)) {

      cat_data_dummy <- recipes::bake(prep$dummy_recipe, new_data = NULL)

      if (is.null(validate_x)) {
        cont_dist_mat <- ndist(
          x = cont_data,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling_cont
        ) |>
          as.matrix()

        cat_dist_mat <- ndist(
          x = cat_data_dummy,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling_cont
        ) |>
          as.matrix()

      } else {
        cont_dist_mat <- ndist(
          x = cont_data,
          validate_x = cont_data_val,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling_cont
        ) |>
          as.matrix()

        cat_data_val_dummy <- recipes::bake(prep$dummy_recipe, new_data = cat_data_val)

        cat_dist_mat <- ndist(
          x = cat_data_dummy,
          validate_x = cat_data_val_dummy,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling_cont
        ) |>
          as.matrix()
      }

      distance_mat <- sqrt((cat_dist_mat^2) + (cont_dist_mat^2))

    } else {
      stop("No variables available to compute distance.", call. = FALSE)
    }

  } else if (prep$preset == "gudmm") {

    no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
    df <- if (!is.null(cont_data) && !is.null(cat_data)) {
      cbind(cont_data, cat_data)
    } else if (!is.null(cont_data)) {
      cont_data
    } else if (!is.null(cat_data)) {
      cat_data
    } else {
      stop("No variables available to compute distance.", call. = FALSE)
    }

    if (is.null(validate_x)) {
      X_matrix <- gudmm_preprocessing(df, no_f_cont)
      Di <- gudmm_distance_dependency_mixed_matrix(X_matrix, no_f_cont, no_f_ord = 0, DM = "DM5")
      distance_mat <- as.matrix(Di)
    } else {
      stop("train to test distances not implemented for this method", call. = FALSE)
    }

  } else if (prep$preset == "dkss") {

    no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
    df <- if (!is.null(cont_data) && !is.null(cat_data)) {
      cbind(cont_data, cat_data)
    } else if (!is.null(cont_data)) {
      cont_data
    } else if (!is.null(cat_data)) {
      cat_data
    } else {
      stop("No variables available to compute distance.", call. = FALSE)
    }

    if (is.null(validate_x)) {
      df <- dkss_preprocessing(df, no_f_cont)
      dkss_result <- kdml::dkss(
        df = df, bw = "np",
        cFUN = "c_gaussian", uFUN = "u_aitken", oFUN = "o_wangvanryzin",
        stan = TRUE, verbose = FALSE
      )
      distance_mat <- dkss_result$distances
    } else {
      stop("train to test distances not implemented for this method", call. = FALSE)
    }

  } else if (prep$preset == "mod_gower") {

    if (is.null(validate_x)) {
      df <- if (!is.null(cont_data) && !is.null(cat_data)) {
        cbind(cont_data, cat_data)
      } else if (!is.null(cont_data)) {
        cont_data
      } else if (!is.null(cat_data)) {
        cat_data
      } else {
        stop("No variables available to compute distance.", call. = FALSE)
      }

      distance_mat <- mg_gower_mod_matrix(df, use_weights = TRUE) |>
        as.matrix()
    } else {
      stop("train to test distances not implemented for this method", call. = FALSE)
    }

  } else {
    stop("Unknown preset: ", prep$preset, call. = FALSE)
  }

  distance_mat
}


#' Mixed-type distance for mixed data
#'
#' Computes a distance or dissimilarity object for mixed-type data,
#' combining continuous and categorical variables.
#'
#' @param x A data frame or matrix of predictors.
#' @param new_data Optional new data to compute distances to \code{x}.
#' @param response Optional response variable used for supervised
#'   categorical distances.
#' @param distance_cont Character string specifying the distance
#'   for continuous variables.
#' @param distance_cat Character string specifying the distance
#'   for categorical variables.
#' @param commensurable Logical; whether to apply commensurability
#'   scaling between variable types.
#' @param scaling_cont Character string specifying scaling for
#'   continuous variables.
#' @param ncomp Optional number of components used in dimensional
#'   reduction. If \code{NULL}, no reduction is applied.
#' @param threshold Optional threshold parameter.
#' @param preset Character string specifying a predefined setup.
#' @param interaction Logical; whether to include interaction-based
#'   distances.
#' @param prop_nn proportion of neighbours to consider when measuring interactions
#' @param score classification metric either "ba" (balanced accuracy) or "logloss"
#' @param decision rule when score is set to ba
#' @return A dissimilarity object.
#' @export
mdist <- function(x, new_data = NULL, response = NULL,
                  distance_cont = "manhattan", distance_cat = "tvd",
                  commensurable = TRUE, scaling_cont = "std",
                  ncomp = NULL, threshold = NULL,
                  preset = "custom", interaction = FALSE,
                  prop_nn = 0.1, score = "ba", decision = "prior_corrected") {

  prep <- .prep_mdist(
    x = x,
    response = {{ response }},
    distance_cont = distance_cont,
    distance_cat = distance_cat,
    commensurable = commensurable,
    scaling_cont = scaling_cont,
    ncomp = ncomp,
    threshold = threshold,
    preset = preset,
    interaction = interaction,
    prop_nn = prop_nn,
    score = score,
    decision = decision
  )

  distance_mat <- .apply_mdist(prep, new_data = new_data)
  distance_mat <- .to_dissimilarity(distance_mat)

  params <- list(
    distance_cont = prep$distance_cont,
    distance_cat  = prep$distance_cat,
    scaling_cont  = prep$scaling_cont,
    commensurable = prep$commensurable,
    ncomp         = prep$ncomp,
    threshold     = prep$threshold,
    rectangular   = (nrow(as.matrix(distance_mat)) != ncol(as.matrix(distance_mat))),
    train_n       = nrow(prep$x_train),
    test_n        = if (!is.null(new_data)) nrow(tibble::as_tibble(new_data)) else NULL,
    cat_p         = if (!is.null(prep$cat_data))  ncol(prep$cat_data)  else 0,
    cont_p        = if (!is.null(prep$cont_data)) ncol(prep$cont_data) else 0,
    response_col  = prep$response_col,
    y             = prep$y,
    preprocessor  = prep
  )

  MDist$new(
    distance = distance_mat,
    preset   = prep$preset,
    data     = prep$x_train,
    params   = params
  )
}
