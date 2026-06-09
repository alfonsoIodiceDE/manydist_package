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
    method_cat,
    commensurable,
    method_num,
    interaction
) {

  preset_defs <- list(
    custom = list(),

    gower = list(
      method_cat    = "matching",
      commensurable = FALSE,
      method_num    = "range",
      interaction   = FALSE
    ),
    euclidean = list(
      method_cat    = "dummy",
      commensurable = FALSE,
      method_num    = "std",
      interaction   = FALSE
    ),
    unbiased_dependent = list(
      method_cat    = "tvd",
      commensurable = TRUE,
      method_num    = "pc_scores",
      interaction   = FALSE
    ),

    hl = list(
      method_cat    = "HLeucl",
      commensurable = FALSE,
      method_num    = "std",
      interaction   = FALSE
    ),

    u_dep = list(
      method_cat    = "tvd",
      commensurable = TRUE,
      method_num    = "pc_scores",
      interaction   = FALSE
    ),

    u_indep = list(
      method_cat    = "matching",
      commensurable = TRUE,
      method_num    = "std",
      interaction   = FALSE
    ),

    u_mix = list(
      method_cat    = "tvd",
      commensurable = TRUE,
      method_num    = "std",
      interaction   = FALSE
    )
  )

  user_args <- list(
    method_cat    = method_cat,
    commensurable = commensurable,
    method_num    = method_num,
    interaction   = interaction
  )

  default_args <- list(
    method_cat    = "tvd",
    commensurable = TRUE,
    method_num    = "std",
    interaction   = FALSE
  )

  if (!preset %in% names(preset_defs)) {
    return(user_args)
  }

  if (identical(preset, "custom")) {
    return(user_args)
  }

  non_default_args <- names(user_args)[
    !vapply(
      names(user_args),
      function(arg) identical(user_args[[arg]], default_args[[arg]]),
      logical(1)
    )
  ]

  if (length(non_default_args) > 0) {
    warning(
      "When `preset` is not 'custom', distance-related arguments are ignored: ",
      paste0("`", non_default_args, "`", collapse = ", "),
      ". Set `preset = 'custom'` to specify them manually.",
      call. = FALSE
    )
  }

  defs <- preset_defs[[preset]]

  list(
    method_cat    = defs$method_cat    %||% method_cat,
    commensurable = defs$commensurable %||% commensurable,
    method_num    = defs$method_num    %||% method_num,
    interaction   = defs$interaction   %||% interaction
  )
}


.normalize_preset <- function(preset) {
  if (preset %in% c("euclidean_onehot", "euclidean_one_hot")) {
    warning(
      "`preset = \"", preset, "\"` is deprecated. ",
      "Use `preset = \"euclidean\"` instead.",
      call. = FALSE
    )
    return("euclidean")
  }

  preset
}

.numeric_metric_from_preset <- function(preset) {
  if (preset %in% c("hl", "euclidean")) {
    "euclidean"
  } else {
    "manhattan"
  }
}

.mdist_generic <- function(
    cont_data, cat_data,
    cont_data_val = NULL, cat_data_val = NULL,
    y = NULL,
    preset, method_cat,
    commensurable, method_num,
    ncomp, threshold,
    interaction, prop_nn, score, decision
) {

  cont_dist_mat <- NULL
  cat_dist_mat  <- NULL

  method <- .numeric_metric_from_preset(preset)

  if (!is.null(cont_data)) {
    cont_dist_mat <- ndist(
      x = cont_data,
      validate_x = cont_data_val,
      method = method,
      commensurable = commensurable,
      scaling = method_num,
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
      method = method_cat,
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
      preset %in% c("hl", "euclidean")) {
    distance_mat <- sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
  }

  distance_mat
}


.prep_mdist <- function(x, response = NULL,
                        method_cat = "tvd",
                        commensurable = TRUE,
                        method_num = "std",
                        ncomp = NULL, threshold = NULL,
                        preset = "custom", interaction = FALSE,
                        prop_nn = 0.1, score = "ba", decision = "prior_corrected",
                        gower_average = TRUE) {

  preset <- .normalize_preset(preset)

  if (identical(method_cat, "tot_var_dist")) {
    method_cat <- "tvd"
  }


  resp_name <- NULL
  response_quo <- rlang::enquo(response)

  if (!rlang::quo_is_null(response_quo)) {
    response_expr <- rlang::quo_get_expr(response_quo)

    if (rlang::is_string(response_expr)) {
      resp_name <- response_expr

      if (!resp_name %in% colnames(x)) {
        stop("`response` must name a column inside `x`.", call. = FALSE)
      }

    } else {
      resp_idx <- tidyselect::eval_select(response_quo, x)

      if (length(resp_idx) != 1L) {
        stop("`response` must select exactly one column.", call. = FALSE)
      }

      resp_name <- names(resp_idx)
    }

    y <- x[[resp_name]]
    x <- x[, setdiff(colnames(x), resp_name), drop = FALSE]
  } else {
    y <- NULL
  }

  x <- tibble::as_tibble(x)

  if (is.null(ncomp)) ncomp <- ncol(x)

  args <- .resolve_preset(
    preset = preset,
    method_cat = method_cat,
    commensurable = commensurable,
    method_num = method_num,
    interaction = interaction
  )

  method_cat    <- args$method_cat
  commensurable <- args$commensurable
  method_num    <- args$method_num
  interaction   <- args$interaction


  cat_data  <- x |> dplyr::select(where(is.factor))
  cont_data <- x |> dplyr::select(where(is.numeric))

  if (ncol(cat_data) == 0) cat_data <- NULL
  if (ncol(cont_data) == 0) cont_data <- NULL

  if (!is.null(cat_data) && identical(method_cat, "tvd") && ncol(cat_data) == 1) {
    warning("'tvd' requires >1 categorical variable. Switching to 'matching'.", call. = FALSE)
    method_cat <- "matching"
  }

  if (!is.null(cont_data) && identical(method_num, "pc_scores") && ncol(cont_data) == 1) {
    warning(
      "With 1 continuous variable, 'pc_scores' is equivalent to standardization. Switching to method_num='std'.",
      call. = FALSE
    )
    method_num <- "std"
  }

  cat_levels <- NULL
  if (!is.null(cat_data)) {
    cat_levels <- purrr::map(cat_data, levels)
  }

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

  if (preset == "euclidean" && !is.null(cat_data)) {
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
      method_cat     = method_cat,
      method_num     = method_num,
      commensurable  = commensurable,
      ncomp          = ncomp,
      threshold      = threshold,
      interaction    = interaction,
      prop_nn        = prop_nn,
      score          = score,
      decision       = decision,
      gower_average  = gower_average,
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

  generic_presets <- c(
    "custom",
    "unbiased_dependent",
    "u_dep",
    "u_indep",
    "u_mix",
    "hl"
  )

  if (prep$preset %in% generic_presets) {

    distance_mat <- .mdist_generic(
      cont_data = cont_data,
      cat_data = cat_data,
      cont_data_val = cont_data_val,
      cat_data_val = cat_data_val,
      y = prep$y,
      preset = prep$preset,
      method_cat = prep$method_cat,
      commensurable = prep$commensurable,
      method_num = prep$method_num,
      ncomp = prep$ncomp,
      threshold = prep$threshold,
      interaction = prep$interaction,
      prop_nn = prep$prop_nn,
      score = prep$score,
      decision = prep$decision
    )

  } else if (prep$preset == "gower") {

    if (is.null(validate_x)) {

      df <- if (!is.null(cont_data) && !is.null(cat_data)) {
        prep$x_train
      } else if (!is.null(cont_data)) {
        cont_data
      } else if (!is.null(cat_data)) {
        cat_data
      } else {
        stop("No variables available to compute distance.", call. = FALSE)
      }

      if (isFALSE(prep$commensurable)) {

        distance_mat <- cluster::daisy(df, metric = "gower") |>
          as.matrix()

        if (isFALSE(prep$gower_average)) {
          distance_mat <- ncol(df) * distance_mat
        }

      } else {

        gowerlist <- df |>
          purrr::map(~ cluster::daisy(tibble::as_tibble(.x), metric = "gower") |>
                       as.matrix())

        gowerlist <- tibble::tibble(gowdist = gowerlist) |>
          dplyr::mutate(commgow = purrr::map(.data$gowdist, ~ .safe_comm(.x)))

        distance_mat <- Reduce(`+`, gowerlist$commgow)

        if (isTRUE(prep$gower_average)) {
          distance_mat <- distance_mat / length(gowerlist$commgow)
        }
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

      if (isTRUE(prep$gower_average)) {
        distance_mat <- distance_mat / length(gowerlist)
      }
    }
  }else if (prep$preset == "euclidean") {

    method <- "euclidean"
    commensurable <- prep$commensurable
    method_num <- prep$method_num

    if (is.null(cat_data) && !is.null(cont_data)) {

      distance_mat <- ndist(
        x = cont_data,
        validate_x = cont_data_val,
        method = method,
        commensurable = commensurable,
        scaling = method_num
      ) |>
        as.matrix()

    } else if (!is.null(cat_data) && is.null(cont_data)) {

      cat_data_dummy <- recipes::bake(prep$dummy_recipe, new_data = NULL)

      if (is.null(validate_x)) {
        distance_mat <- ndist(
          x = cat_data_dummy,
          method = method,
          commensurable = commensurable,
          scaling = method_num
        ) |>
          as.matrix()
      } else {
        cat_data_val_dummy <- recipes::bake(prep$dummy_recipe, new_data = cat_data_val)

        distance_mat <- ndist(
          x = cat_data_dummy,
          validate_x = cat_data_val_dummy,
          method = method,
          commensurable = commensurable,
          scaling = method_num
        ) |>
          as.matrix()
      }

    } else if (!is.null(cat_data) && !is.null(cont_data)) {

      cat_data_dummy <- recipes::bake(prep$dummy_recipe, new_data = NULL)

      if (is.null(validate_x)) {
        cont_dist_mat <- ndist(
          x = cont_data,
          method = method,
          commensurable = commensurable,
          scaling = method_num
        ) |>
          as.matrix()

        cat_dist_mat <- ndist(
          x = cat_data_dummy,
          method = method,
          commensurable = commensurable,
          scaling = method_num
        ) |>
          as.matrix()

      } else {
        cont_dist_mat <- ndist(
          x = cont_data,
          validate_x = cont_data_val,
          method = method,
          commensurable = commensurable,
          scaling = method_num
        ) |>
          as.matrix()

        cat_data_val_dummy <- recipes::bake(prep$dummy_recipe, new_data = cat_data_val)

        cat_dist_mat <- ndist(
          x = cat_data_dummy,
          validate_x = cat_data_val_dummy,
          method = method,
          commensurable = commensurable,
          scaling = method_num
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


#' Mixed-type dissimilarities for distance-based learning
#'
#' Computes a dissimilarity object for numerical, categorical, or mixed-type
#' data. The function combines continuous and categorical components according
#' to either a predefined `preset` or a user-defined custom specification.
#'
#' `mdist()` is the main distance-construction function in `manydist`. It can
#' return ordinary train-train dissimilarities or rectangular test-to-training
#' dissimilarities when `new_data` is supplied. The resulting object stores both
#' the dissimilarity matrix and metadata about the distance specification that
#' was used.
#'
#' @param x A data frame or matrix containing the training observations. Columns
#'   can be numeric, factors, or a mixture of both.
#' @param new_data Optional data frame or matrix containing new observations.
#'   If supplied, distances are computed from rows of `new_data` to rows of
#'   `x`, producing a rectangular test-to-training dissimilarity matrix.
#' @param response Optional response variable used for response-aware
#'   categorical dissimilarities. It can be supplied as an unquoted column name
#'   or as a character string. The response column is removed from the
#'   predictors before computing distances.
#' @param method_cat Character string specifying the categorical-variable
#'   dissimilarity used when `preset = "custom"`. Common values include
#'   `"matching"` and `"tvd"`. Use [all_dist_method_specs()] to inspect
#'   available methods.
#' @param method_num Character string specifying the numerical-variable
#'   preprocessing used when `preset = "custom"`. Available options include
#'   `"none"` for no preprocessing, `"std"` for standard-deviation scaling,
#'   `"range"` for range scaling, `"robust"` for inter-quartile-range-based
#'   scaling, and `"pc_scores"` for principal-component score scaling.
#' @param commensurable Logical. If `TRUE`, dissimilarities are scaled so that
#'   the average contribution of each variable to the overall distance is equal
#'   to 1.
#' @param ncomp Integer or `NULL`. Number of principal components to retain
#'   when `method_num = "pc_scores"`. If `NULL`, all available components are
#'   used unless `threshold` is supplied and supported by the underlying method.
#' @param threshold Numeric or `NULL`. Optional cumulative variance threshold
#'   used when `method_num = "pc_scores"`.
#' @param preset Character string specifying a predefined distance
#'   specification. Available values include `"custom"`, `"gower"`,
#'   `"unbiased_dependent"`, `"u_dep"`, `"u_indep"`, `"u_mix"`, `"hl"`,
#'   `"gudmm"`, `"dkss"`, `"mod_gower"`, and `"euclidean"`.
#'   When `preset` is not `"custom"`, arguments such as `method_cat`,
#'   `method_num`, `commensurable`, and `interaction` are handled by the preset
#'   and user-supplied values for those arguments are ignored.
#' @param interaction Logical. If `TRUE`, adds an interaction-aware
#'   continuous-categorical component based on local predictive separability.
#' @param prop_nn Numeric. Proportion of nearest neighbours used when
#'   `interaction = TRUE`.
#' @param score Character string specifying the score used when
#'   `interaction = TRUE`. Available values include `"ba"` for balanced
#'   accuracy and `"logloss"`.
#' @param decision Character string specifying the decision rule used when
#'   `score = "ba"`. The default is `"prior_corrected"`.
#' @param gower_average Logical; only used when `preset = "gower"`. If `TRUE`,
#'   returns the standard Gower dissimilarity averaged over variables, matching
#'   the scale of [cluster::daisy()] with `metric = "gower"`. If `FALSE`,
#'   returns the sum of per-variable Gower contributions, equivalent to
#'   multiplying the averaged Gower dissimilarity by the number of active
#'   variables.
#'
#' @details
#' With `preset = "custom"`, users manually choose the numerical preprocessing,
#' categorical dissimilarity, commensurability, and optional interaction term.
#'
#' The `"gower"` preset follows the usual Gower construction based on range
#' scaling for continuous variables and matching dissimilarities for categorical
#' variables. The `gower_average` argument controls whether the result is
#' averaged over variables or returned as a sum of variable-wise contributions.
#'
#' The `"u_dep"`, `"unbiased_dependent"`, `"u_indep"`, and `"u_mix"` presets are
#' convenience specifications for unbiased or commensurable mixed-variable
#' dissimilarities. The `"euclidean"` preset computes a Euclidean
#' distance after one-hot encoding categorical variables. The `"gudmm"`,
#' `"dkss"`, and `"mod_gower"` presets provide additional mixed-type distance
#' constructions. Some presets currently support only train-train distances and
#' will stop if `new_data` is supplied.
#'
#' Use [all_dist_method_specs()] to inspect the available distance components
#' and method specifications.
#'
#' @return An object of class `"MDist"`. The object contains the computed
#'   dissimilarity in its `$distance` field, the selected `preset`, the training
#'   data, and a list of parameters describing the fitted distance
#'   specification. Square train-train dissimilarities are stored as
#'   `"dissimilarity"`/`"dist"` objects; rectangular test-to-training
#'   dissimilarities are stored as `"dissimilarity"`/`"matrix"` objects.
#'
#' @seealso [step_mdist()], [all_dist_method_specs()]
#'
#' @examples
#' if (requireNamespace("palmerpenguins", quietly = TRUE)) {
#'   data("penguins", package = "palmerpenguins")
#'
#'   penguins_small <- palmerpenguins::penguins |>
#'     dplyr::select(
#'       bill_length_mm, bill_depth_mm, flipper_length_mm,
#'       body_mass_g, species, island, sex
#'     ) |>
#'     tidyr::drop_na()
#'
#'   # Gower distance on mixed-type data
#'   d_gower <- mdist(penguins_small, preset = "gower")
#'   d_gower
#'
#'   # Custom mixed-type specification
#'   d_custom <- mdist(
#'     penguins_small,
#'     preset = "custom",
#'     method_cat = "matching",
#'     method_num = "std",
#'     commensurable = TRUE
#'   )
#'
#'   d_custom
#'
#'   # Train-to-new-data distances
#'   penguin_split <- rsample::initial_split(penguins_small, prop = 0.75)
#'   penguin_train <- rsample::training(penguin_split)
#'   penguin_test  <- rsample::testing(penguin_split)
#'
#'   d_new <- mdist(
#'     penguin_train,
#'     new_data = penguin_test,
#'     preset = "gower"
#'   )
#'
#'   d_new
#' }
#'
#' @export
mdist <- function(x, new_data = NULL, response = NULL,
                  method_cat = "tvd", method_num = "std",
                  commensurable = TRUE,
                  ncomp = NULL, threshold = NULL,
                  preset = "custom", interaction = FALSE,
                  prop_nn = 0.1, score = "ba", decision = "prior_corrected",
                  gower_average = TRUE) {

  prep <- .prep_mdist(
    x = x,
    response = {{ response }},
    method_cat = method_cat,
    commensurable = commensurable,
    method_num = method_num,
    ncomp = ncomp,
    threshold = threshold,
    preset = preset,
    interaction = interaction,
    prop_nn = prop_nn,
    score = score,
    decision = decision,
    gower_average = gower_average
  )

  distance_mat <- .apply_mdist(prep, new_data = new_data)
  distance_mat <- .to_dissimilarity(distance_mat)

  params <- list(
    method_cat  = prep$method_cat,
    method_num  = prep$method_num,
    commensurable = prep$commensurable,
    gower_average = prep$gower_average,
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
