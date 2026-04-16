cdist <- function(x, response = NULL, validate_x = NULL,
                  method = "tvd", commensurable = TRUE, weights = 1) {

  full_delta  <- NULL
  level_stop  <- NULL
  level_start <- NULL
  delta_tmp   <- NULL
  all_of      <- NULL
  map2        <- NULL
  .x          <- NULL
  a           <- NULL
  b           <- NULL
  blocks      <- NULL
  id          <- NULL
  delta_ms    <- NULL

  # ------------------------------------------------------------------
  # helpers
  # ------------------------------------------------------------------

  .dummy_encode <- function(train_x, new_x = NULL, nominal_only = FALSE) {
    rec <- recipes::recipe(~ ., data = tibble::as_tibble(train_x))

    if (isTRUE(nominal_only)) {
      rec <- rec |>
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = TRUE)
    } else {
      rec <- rec |>
        recipes::step_dummy(recipes::all_predictors(), one_hot = TRUE)
    }

    rec <- recipes::prep(rec, training = tibble::as_tibble(train_x))

    if (is.null(new_x)) {
      recipes::bake(rec, new_data = NULL) |> as.matrix()
    } else {
      recipes::bake(rec, new_data = tibble::as_tibble(new_x)) |> as.matrix()
    }
  }

  .expand_weights <- function(weights, Q) {
    weightsexp <- NULL
    for (i in seq_along(Q)) {
      weightsexp <- c(weightsexp, rep(weights[i], Q[i]))
    }
    diag(weightsexp, nrow = length(weightsexp), ncol = length(weightsexp))
  }

  .assemble_dist <- function(Z_train, delta, validate_x = NULL, x_train = NULL,
                             weights = 1, Q = NULL, nominal_only = FALSE) {
    if (is.null(validate_x)) {
      Z_val <- Z_train
    } else {
      Z_val <- .dummy_encode(train_x = x_train, new_x = validate_x, nominal_only = nominal_only)
    }

    if (length(weights) == 1 && is.null(dim(weights))) {
      return(Z_val %*% delta %*% t(Z_train))
    }

    if (!is.null(dim(weights))) {
      weights <- diag(weights)
    }

    W <- .expand_weights(weights = weights, Q = Q)
    Z_val %*% W %*% delta %*% W %*% t(Z_train)
  }

  # ------------------------------------------------------------------
  # normalize inputs
  # ------------------------------------------------------------------

  if (length(method) == 1 && method %in% c("tot_var_dist", "tvd")) {
    method <- "tvd"
  } else if (length(method) > 1) {
    method[method %in% c("tot_var_dist", "tvd")] <- "tvd"
  }

  y <- response
  x <- tibble::as_tibble(x) |>
    dplyr::mutate(dplyr::across(where(is.factor), forcats::fct_drop))

  response_irrelevant_methods <- c(
    "none", "st_dev", "HL", "cat_dis", "HLeucl", "mca",
    "matching", "eskin", "goodall_3", "goodall_4",
    "iof", "of", "lin", "var_entropy", "var_mutability"
  )

  profile_based_methods <- c(
    "tvd", "le_and_ho",
    philentropy::getDistMethods()
  )

  multi_method_not_supported <- c(
    "none", "st_dev", "HL", "cat_dis", "HLeucl", "mca"
  )

  if (!is.null(y)) {
    ignored_response <- unique(method[method %in% response_irrelevant_methods])

    if (length(ignored_response) > 0) {
      warning(
        "For method(s) ",
        paste(sprintf("'%s'", ignored_response), collapse = ", "),
        ", category dissimilarities do not depend on conditional profiles; ",
        "`response` was therefore ignored.",
        call. = FALSE
      )
    }
  }

  if ((length(method) > 1) && any(method %in% multi_method_not_supported)) {
    stop(
      "Specifying a different method for each variable is not currently possible for 'st_dev', 'HL', 'cat_dis', 'HLeucl', or 'mca'.",
      call. = FALSE
    )
  }
  # ------------------------------------------------------------------
  # single method
  # ------------------------------------------------------------------

  if (length(method) == 1) {

    if (method %in% c("none", "st_dev", "HL", "cat_dis", "HLeucl", "mca")) {

      outindbased <- indicator_based(
        x = x,
        validate_x = validate_x,
        commensurable = commensurable,
        scaling = method,
        weights = 1
      )

      distance_mat <- outindbased$distance_mat
      delta        <- outindbased$delta
      delta_ms     <- outindbased$delta_ms
      delta_names  <- outindbased$delta_names

    } else {

      out_delta   <- cat_delta(x = x, y = y, method_cat = method)
      delta       <- out_delta[[method]] |> data.matrix()
      delta_names <- out_delta$delta_names
      Z           <- out_delta$Z |> data.matrix()

      colnames(delta) <- delta_names
      rownames(delta) <- delta_names

      if (is.null(dim(x))) {
        Q <- nlevels(x)
      } else {
        Q <- purrr::map_dbl(x, nlevels)
      }

      if (isTRUE(commensurable)) {

        distance_mat <- commensurability_for_cat(
          x = x, y = y, validate_x = validate_x, delta = delta
        )

      } else {

        distance_mat <- .assemble_dist(
          Z_train     = Z,
          delta       = delta,
          validate_x  = validate_x,
          x_train     = x,
          weights     = weights,
          Q           = Q,
          nominal_only = is.null(validate_x)
        )
      }
    }

  } else {

    # ----------------------------------------------------------------
    # different method per variable
    # ----------------------------------------------------------------

    method_vec <- method

    if (is.null(dim(x))) {
      Q <- nlevels(x)
    } else {
      Q <- purrr::map_dbl(x, nlevels)
    }

    nvar      <- length(Q)
    level_pos <- data.table::data.table(
      start = c(1, cumsum(Q)[-length(Q)] + 1),
      stop  = cumsum(Q)
    )

    delta_structure <- tibble::tibble(method = method_vec) |>
      dplyr::mutate(
        x = purrr::map(.x = method_vec, ~ tibble::as_tibble(x)),
        delta_tmp = purrr::map2(.x = x, .y = method, .f = ~ cat_delta(x = .x, y = y, method_cat = .y)),
        delta_names = purrr::map(.x = delta_tmp, .f = ~ .x$delta_names),
        full_delta = purrr::map2(.x = delta_tmp, .y = method, .f = ~ .x[[.y]]),
        level_start = level_pos$start,
        level_stop  = level_pos$stop,
        delta_block = purrr::pmap(
          .l = list(a = full_delta, b = level_start, c = level_stop),
          .f = function(a, b, c) a[b:c, b:c]
        )
      )

    Z           <- delta_structure$delta_tmp[[1]]$Z |> as.matrix()
    delta       <- Matrix::bdiag(delta_structure$delta_block) |> as.matrix()
    delta_names <- delta_structure$delta_names |> unlist()

    if (isTRUE(commensurable)) {

      distance_mat <- commensurability_for_cat(
        x = x, y = y, validate_x = validate_x, delta = delta
      )

    } else {

      distance_mat <- .assemble_dist(
        Z_train     = Z,
        delta       = delta,
        validate_x  = validate_x,
        x_train     = x,
        weights     = weights,
        Q           = Q,
        nominal_only = FALSE
      )
    }
  }

  # ------------------------------------------------------------------
  # output
  # ------------------------------------------------------------------

  out_catdist <- list()
  out_catdist$distance_mat <- .to_dissimilarity(distance_mat)
  out_catdist$delta        <- delta
  out_catdist$delta_ms     <- delta_ms
  out_catdist$delta_names  <- delta_names

  out_catdist
}
