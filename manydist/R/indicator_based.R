indicator_based <- function(x, validate_x, commensurable = FALSE, scaling = "none", weights = 1) {
  factor_name = NULL
  delta = NULL
  by_var_dist = NULL
  mean_by_var_dist = NULL
  .x = NULL
  eta = NULL
  scaled_Zs = NULL
  comm_dist = NULL
  scaled_val_Zs = NULL

  delta_output = cat_delta(x = x, method_cat = "matching")
  Z = delta_output$Z
  delta_m = delta_output$matching
  delta_ms <- NULL
  delta_names = delta_output$delta_names
  p <- ncol(x)

  if (scaling == "none") {
    delta_m = 2 * delta_m
    if (!commensurable) {
      if (is.null(validate_x)) {
        cat_dist_mat = Z %*% delta_m %*% t(Z)
      } else {
        validate_Z = as.matrix(bake(prep(step_dummy(recipe(validate_x, ~.),
                                                    all_predictors(), one_hot = TRUE), training = as_tibble(x)),
                                    new_data = validate_x))
        cat_dist_mat = validate_Z %*% delta_m %*% t(Z)
      }
    } else {
      prep_list = map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                          all_predictors(), one_hot = TRUE), training = as_tibble(.x)))
      Z_list = map(prep_list, ~bake(.x, new_data = NULL))
      Q = map_dbl(x, nlevels)
      levels_identifier = rep(names(Q), times = as.vector(Q))
      if (is.null(validate_x)) {
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_m[levels_identifier == .x,
                                                                                               levels_identifier == .x])),
                                              Zs = Z_list,
                                              by_var_dist = map2(.x = Z_list, .y = delta,
                                                                 ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      } else {
        validate_Z_list = map2(.x = validate_x, .y = prep_list,
                               ~bake(.y, new_data = as_tibble(.x)))
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_m[levels_identifier == .x,
                                                                                               levels_identifier == .x])),
                                              Zs = Z_list,
                                              val_Zs = validate_Z_list,
                                              by_var_dist = pmap(.l = list(a = Z_list, b = delta, c = validate_Z_list),
                                                                 .f = function(a, b, c) {
                                                                   return(as.matrix(c) %*% b %*% t(as.matrix(a)))
                                                                 }),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      }
      cat_dist_mat = Reduce(`+`, pull(commensurable_dist_structure, comm_dist))
    }
  } else if (scaling == "st_dev") {
    qs_vec = unlist(map(as_tibble(x), function(x) {
      as.vector(rep(nlevels(x), nlevels(x)))
    }))
    qs_diag = diag(1/sqrt(unlist(qs_vec)))
    prep_c_Z = prep(step_center(step_dummy(recipe(as_tibble(x), ~.),
                                           all_predictors(), one_hot = TRUE), all_predictors()), training = x)
    c_Z = as.matrix(bake(prep_c_Z, new_data = NULL))
    S = (1/nrow(c_Z)) * t(c_Z) %*% c_Z
    inv_sq_S_d = diag(1/sqrt(diag(S)))
    prep_Z = prep(step_dummy(recipe(x, ~.), all_nominal(), one_hot = TRUE), training = x)
    Z = as.matrix(bake(prep_Z, new_data = NULL))
    Z = apply(Z, 2, as.numeric)
    delta_ms <- inv_sq_S_d %*% qs_diag %*% delta_m %*% inv_sq_S_d
    if (!commensurable) {
      if (is.null(validate_x)) {
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
      } else {
        val_Z = as.matrix(bake(prep_Z, new_data = validate_x))
        val_Z = apply(val_Z, 2, as.numeric)
        cat_dist_mat = val_Z %*% delta_ms %*% t(Z)
      }
    } else {
      Z_prep_list = map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                            all_predictors(), one_hot = TRUE), training = as_tibble(.x)))
      Z_list = map(Z_prep_list, ~bake(.x, new_data = NULL))
      Q = map_dbl(x, nlevels)
      levels_identifier = rep(names(Q), times = as.vector(Q))
      if (is.null(validate_x)) {
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_ms[levels_identifier == .x,
                                                                                                levels_identifier == .x])),
                                              Zs = Z_list,
                                              by_var_dist = map2(.x = Z_list, .y = delta,
                                                                 ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
        cat_dist_mat = Reduce(`+`, pull(commensurable_dist_structure, comm_dist))
      } else {
        val_Z_list = map2(.x = Z_prep_list, .y = validate_x,
                          ~bake(.x, new_data = as_tibble(.y)))
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_ms[levels_identifier == .x,
                                                                                                levels_identifier == .x])),
                                              Zs = Z_list,
                                              val_Zs = val_Z_list,
                                              by_var_dist = pmap(.l = list(a = Z_list, b = delta, c = val_Z_list),
                                                                 function(a, b, c) as.matrix(c) %*% b %*% t(as.matrix(a))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
        cat_dist_mat = Reduce(`+`, pull(commensurable_dist_structure, comm_dist))
      }
    }
  } else if (scaling == "HL") {
    eta_vec = map(as_tibble(x), function(x = .x) {
      as.vector(rep(fpc::distancefactor(cat = nlevels(x), catsizes = table(x)), nlevels(x)))
    })
    eta_diag = diag(unlist(eta_vec))
    delta_ms <- 2 * delta_m %*% eta_diag

    # FIX: Create prep_Z once at the start of this branch
    prep_Z = prep(step_dummy(recipe(x, ~.), all_nominal(), one_hot = TRUE), training = x)
    Z = as.matrix(bake(prep_Z, new_data = NULL))
    Z = apply(Z, 2, as.numeric)

    if (!commensurable) {
      if (is.null(validate_x)) {
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
      } else {
        val_Z = as.matrix(bake(prep_Z, new_data = validate_x))
        val_Z = apply(val_Z, 2, as.numeric)
        cat_dist_mat = val_Z %*% delta_ms %*% t(Z)
      }
    } else {
      if (is.null(validate_x)) {
        prep_Z_list = map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                              all_predictors(), one_hot = TRUE), training = as_tibble(.x)))
        Z_list = map(prep_Z_list, ~bake(.x, new_data = NULL))
        Q = map_dbl(x, nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_ms[levels_identifier == .x,
                                                                                                levels_identifier == .x])),
                                              by_var_dist = map2(.x = Z_list, .y = delta,
                                                                 ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      } else {
        prep_Z_list = map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                              all_predictors(), one_hot = TRUE), training = as_tibble(.x)))
        Z_list = map(prep_Z_list, ~bake(.x, new_data = NULL))
        val_Z_list = map2(.x = prep_Z_list, .y = validate_x,
                          ~bake(.x, new_data = as_tibble(.y)))
        Q = map_dbl(x, nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_ms[levels_identifier == .x,
                                                                                                levels_identifier == .x])),
                                              val_Zs = val_Z_list,
                                              by_var_dist = pmap(.l = list(a = Z_list, b = delta, c = val_Z_list),
                                                                 function(a, b, c) as.matrix(c) %*% b %*% t(as.matrix(a))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      }
      cat_dist_mat = Reduce(`+`, pull(commensurable_dist_structure, comm_dist))
    }
  } else if (scaling == "cat_dis") {
    qs_vec = unlist(map(as_tibble(x), function(x) {
      as.vector(rep(nlevels(x), nlevels(x)))
    }))
    qs_diag = diag(1/unlist(qs_vec))
    c_Z = as.matrix(bake(prep(step_center(step_dummy(recipe(as_tibble(x), ~.),
                                                     all_predictors(), one_hot = TRUE), all_predictors()), training = x), new_data = NULL))
    S = (1/nrow(c_Z)) * t(c_Z) %*% c_Z
    inv_sq_S_d = diag(1/sqrt(diag(S)))
    prep_Z = prep(step_dummy(recipe(x, ~.), all_nominal(), one_hot = TRUE), training = x)
    Z = as.matrix(bake(prep_Z, new_data = NULL))
    Z = apply(Z, 2, as.numeric)
    Zs = Z %*% inv_sq_S_d %*% diag(1/sqrt(unlist(qs_vec)))
    delta_ms = (qs_diag %*% inv_sq_S_d %*% delta_m %*% inv_sq_S_d)
    if (!commensurable) {
      if (is.null(validate_x)) {
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
      } else {
        val_Z = as.matrix(bake(prep_Z, new_data = validate_x))
        val_Z = apply(val_Z, 2, as.numeric)
        cat_dist_mat = val_Z %*% delta_ms %*% t(Z)
      }
    } else {
      if (is.null(validate_x)) {
        prep_Z_list = map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                              all_predictors(), one_hot = TRUE), training = as_tibble(.x)))
        Z_list = map(prep_Z_list, ~bake(.x, new_data = NULL))
        Q = map_dbl(x, nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_ms[levels_identifier == .x,
                                                                                                levels_identifier == .x])),
                                              by_var_dist = map2(.x = Z_list, .y = delta,
                                                                 ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      } else {
        prep_Z_list = map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                              all_predictors(), one_hot = TRUE), training = as_tibble(.x)))
        Z_list = map(prep_Z_list, ~bake(.x, new_data = NULL))
        val_Z_list = map2(.x = prep_Z_list, .y = validate_x,
                          ~bake(.x, new_data = as_tibble(.y)))
        Q = map_dbl(x, nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_ms[levels_identifier == .x,
                                                                                                levels_identifier == .x])),
                                              val_Zs = val_Z_list,
                                              by_var_dist = pmap(.l = list(a = Z_list, b = delta, c = val_Z_list),
                                                                 function(a, b, c) as.matrix(c) %*% b %*% t(as.matrix(a))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      }
      cat_dist_mat = Reduce(`+`, pull(commensurable_dist_structure, comm_dist))
    }
  } else if (scaling == "HLeucl") {
    delta_ms <- delta_m

    # FIX: Add training = x to prep()
    prep_Z <- prep(step_dummy(step_rm(recipe(x, ~.), all_numeric_predictors()),
                              all_nominal(), one_hot = TRUE), training = x)
    Z <- as.matrix(bake(prep_Z, new_data = NULL))

    eta_vec = map(as_tibble(x), function(x = .x) {
      as.vector(rep(fpc::distancefactor(cat = nlevels(x), catsizes = table(x)), nlevels(x)))
    })
    Zs <- (Z %*% diag(unlist(eta_vec)))

    if (!commensurable) {
      if (is.null(validate_x)) {
        cat_dist_mat = as.matrix(daisy(Zs, metric = "euclidean", warnType = FALSE))
      } else {
        # FIX: Use the already-created prep_Z instead of recreating it
        val_Z = as.matrix(bake(prep_Z, new_data = validate_x))
        val_Zs = (val_Z %*% diag(unlist(eta_vec)))
        cat_dist_mat = as.matrix(Rfast::dista(xnew = val_Zs, x = Zs, type = "euclidean"))
      }
    } else {
      if (is.null(validate_x)) {
        Z_list <- map(x, ~bake(prep(step_dummy(recipe(as_tibble(.x), ~.),
                                               all_nominal(), one_hot = TRUE), training = as_tibble(.x)), new_data = NULL))
        eta_vec <- map(as_tibble(x), function(x) {
          as.vector(rep(fpc::distancefactor(cat = nlevels(x), catsizes = table(x)), nlevels(x)))
        })
        commensurable_dist_structure <- mutate(tibble(factor_name = names(x)),
                                               Zs = Z_list,
                                               eta = eta_vec,
                                               scaled_Zs = map2(Zs, eta, ~as.matrix(.x) %*% diag(.y)),
                                               by_var_dist = map(scaled_Zs, ~as.matrix(daisy(.x, metric = "euclidean", warnType = FALSE))),
                                               mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                               comm_dist = map2(by_var_dist, mean_by_var_dist, ~.x/.y))
      } else {
        prep_Z_list <- map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                               all_nominal(), one_hot = TRUE), training = as_tibble(.x)))
        Z_list = map(prep_Z_list, ~bake(.x, new_data = NULL))
        val_Z_list = map2(.x = prep_Z_list, .y = validate_x,
                          ~bake(.x, new_data = as_tibble(.y)))
        eta_vec <- map(as_tibble(x), function(x) {
          as.vector(rep(fpc::distancefactor(cat = nlevels(x), catsizes = table(x)), nlevels(x)))
        })
        commensurable_dist_structure <- mutate(tibble(factor_name = names(x)),
                                               Zs = Z_list,
                                               val_Zs = val_Z_list,
                                               eta = eta_vec,
                                               scaled_Zs = map2(Zs, eta, ~as.matrix(.x) %*% diag(.y)),
                                               scaled_val_Zs = map2(val_Zs, eta, ~as.matrix(.x) %*% diag(.y)),
                                               by_var_dist = map2(.x = scaled_val_Zs, .y = scaled_Zs,
                                                                  ~as.matrix(Rfast::dista(xnew = .x, x = .y, type = "euclidean"))),
                                               mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                               comm_dist = map2(by_var_dist, mean_by_var_dist, ~.x/.y))
      }
      cat_dist_mat <- Reduce(`+`, pull(commensurable_dist_structure, comm_dist))
    }
  } else if (scaling == "mca") {
    qs_vec <- unlist(map(as_tibble(dplyr::select(x, where(is.factor))), function(x) {
      as.vector(rep(nlevels(x), nlevels(x)))
    }))
    prep_Z <- prep(step_dummy(step_rm(recipe(x, ~.), all_numeric_predictors()),
                              all_nominal(), one_hot = TRUE), training = x)
    Z = as.matrix(bake(prep_Z, new_data = NULL))
    Z <- apply(Z, 2, as.numeric)
    D_p <- diag(colSums(Z)/nrow(Z))
    inv_sq_D_p <- diag(1/sqrt(diag(D_p)))
    delta_ms <- inv_sq_D_p %*% delta_m %*% inv_sq_D_p
    if (!commensurable) {
      if (is.null(validate_x)) {
        cat_dist_mat = Z %*% inv_sq_D_p %*% delta_m %*% inv_sq_D_p %*% t(Z)
      } else {
        val_Z = as.matrix(bake(prep_Z, new_data = validate_x))
        val_Z <- apply(val_Z, 2, as.numeric)
        cat_dist_mat = val_Z %*% inv_sq_D_p %*% delta_m %*% inv_sq_D_p %*% t(Z)
      }
    } else {
      prep_Z_list = map(x, ~prep(step_dummy(recipe(as_tibble(.x), ~.),
                                            all_predictors(), one_hot = TRUE), training = as_tibble(.x)))
      Z_list = map(prep_Z_list, ~bake(.x, new_data = NULL))
      Q = map_dbl(x, nlevels)
      levels_identifier = rep(names(Q), times = as.vector(Q))
      if (is.null(validate_x)) {
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_m[levels_identifier == .x,
                                                                                               levels_identifier == .x])),
                                              Zs = Z_list,
                                              by_var_dist = map2(.x = Z_list, .y = delta,
                                                                 ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      } else {
        val_Z_list = map2(.x = prep_Z_list, .y = validate_x,
                          ~bake(.x, new_data = as_tibble(.y)))
        commensurable_dist_structure = mutate(tibble(factor_name = names(Q)),
                                              delta = map(.x = factor_name, ~as.matrix(delta_m[levels_identifier == .x,
                                                                                               levels_identifier == .x])),
                                              Zs = Z_list,
                                              val_Zs = val_Z_list,
                                              by_var_dist = pmap(.l = list(a = Z_list, b = delta, c = val_Z_list),
                                                                 .f = function(a, b, c) {
                                                                   return(as.matrix(c) %*% b %*% t(as.matrix(a)))
                                                                 }),
                                              mean_by_var_dist = map_dbl(by_var_dist, ~mean(.x)),
                                              comm_dist = map2(.x = by_var_dist, .y = mean_by_var_dist, ~.x/.y))
      }
      cat_dist_mat = Reduce(`+`, pull(commensurable_dist_structure, comm_dist))
    }
  }

  out_catdist = list()
  out_catdist$distance_mat = cat_dist_mat
  out_catdist$delta = delta_m
  out_catdist$delta_ms = delta_ms
  out_catdist$delta_names = delta_names
  return(out_catdist)
}
