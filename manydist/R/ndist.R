ndist = function (x, validate_x = NULL, commensurable = FALSE, method = "manhattan",
                  sig = NULL, scaling = "none", ncomp = ncol(x), threshold = NULL,
                  weights = rep(1, ncol(x)))
{


  # print(method)
  # print(ncol(x))
  .x <- NULL
  y <- NULL
  by_var_dist_w <- NULL
  rec_x = recipe(~., data = x)
  tr_x = x
  if (scaling == "std") {
    # print("here?")
    rec_x = rec_x |>
      step_scale(all_predictors())

    prepped_rec_x = rec_x |> prep(training=tr_x)

    x = bake(prepped_rec_x, new_data = NULL)

    if (!is.null(validate_x)) {
      validate_x = bake(prepped_rec_x, new_data= validate_x)
    }
  }
  else if (scaling == "range") {

    rec_x = rec_x |> step_range(all_predictors())
    prepped_rec_x = rec_x |> prep(training=tr_x)
    x = bake(prepped_rec_x, new_data = NULL)

    if (!is.null(validate_x)) {
      validate_x = bake(prepped_rec_x,  new_data = validate_x)
    }
  }
  else if (scaling == "pc_scores") {
    if (is.null(threshold)) {
      #  print("ncomp is fixed to ncol(x), has to be changed")
      # ncomp = ncol(x)

      rec_x = rec_x |>
        step_normalize(all_predictors()) |>
        step_pca(all_predictors(), num_comp = ncomp)

      prepped_rec_x = rec_x |> prep(training = tr_x)


      x = bake(prepped_rec_x, new_data = NULL)


      if (!is.null(validate_x)) {
        validate_x = bake(prepped_rec_x, new_data = validate_x)
      }
    }
    else {
      rec_x = rec_x |>
        step_normalize(all_predictors()) |>
        step_pca(all_predictors(), threshold = threshold)

      prepped_rec_x = rec_x |> prep(training = tr_x)

      x = bake(prepped_rec_x, new_data = NULL)

      if (!is.null(validate_x)) {
        validate_x = bake(prepped_rec_x, new_data = validate_x)
      }
    }
  }else if (scaling == "robust") {
    rec_x = rec_x |>
      step_mutate(across(everything(), ~(. - median(.))/IQR(.)))

    prepped_rec_x = rec_x |> prep(training = tr_x)

    x = bake(prepped_rec_x, new_data = NULL)

    if (!is.null(validate_x)) {
      validate_x = bake(prepped_rec_x,
                        new_data = validate_x)
    }
  }
  else if (scaling == "none") {
  }
  if (method == "mahalanobis") {
    distance = mahalanobis_like(x = x, validate_x = validate_x,
                                sig = sig)
  } else {
    if (commensurable == TRUE) {
      if (is.null(validate_x)) {
        by_var_dist = map(.x = as_tibble(x), .f = function(x = .x) {
          b_v_d = daisy(data.frame(x), metric = method,
                        warnBin = FALSE) %>% as.matrix()
          b_v_d = b_v_d/mean(b_v_d)
          return(b_v_d)
        })
      }
      else {
        # save(file = "try_dat.RData",x,validate_x)
        by_var_dist = map2(
          .x = as_tibble(x),
          .y = as_tibble(validate_x),
          .f = function(x, xnew) {
            x    <- as.numeric(x)
            xnew <- as.numeric(xnew)

            # 1D manhattan (and 1D euclidean) distance:
            b_v_d <- abs(outer(xnew, x, "-"))

            b_v_d / mean(b_v_d)
          }
        )

      }
    }
    else if (commensurable == FALSE) {
      if (is.null(validate_x)) {
        by_var_dist = map(.x = as_tibble(x), .f = function(x = .x) {
          b_v_d = daisy(data.frame(x), metric = method,
                        warnBin = FALSE) |> as.matrix()
          return(b_v_d)
        })
      }else {
        by_var_dist = map2(
          .x = as_tibble(x),
          .y = as_tibble(validate_x),
          .f = function(x, xnew) {
            x    <- as.numeric(x)
            xnew <- as.numeric(xnew)

            # 1D manhattan (and 1D euclidean) distance:
            b_v_d <- abs(outer(xnew, x, "-"))

          }
        )


      }
    }
    if (method == "euclidean") {
      by_var_structure = mutate(tibble(by_var_dist = by_var_dist,
                                       weights = weights), by_var_dist_w = map2(.x = by_var_dist,
                                                                                .y = weights, ~(.x^2) * .y))
      distance = as.matrix(Reduce(`+`, pull(by_var_structure,
                                            by_var_dist_w)))
      distance = sqrt(distance)
    }
    else {

      by_var_structure = mutate(tibble(by_var_dist = by_var_dist,
                                       weights = weights), by_var_dist_w = map2(.x = by_var_dist,
                                                                                .y = weights, ~.x * .y))
      distance = as.matrix(Reduce(`+`, pull(by_var_structure,
                                            by_var_dist_w)))
    }
  }

  to_dissimilarity <- function(dist_matrix, reference = NULL) {
    # If square and symmetric → convert to "dist"
    if (nrow(dist_matrix) == ncol(dist_matrix) &&
        isTRUE(all.equal(dist_matrix, t(dist_matrix), tolerance = 1e-10))) {

      d <- as.dist(dist_matrix)
      class(d) <- c("dissimilarity", "dist")

      # Copy optional attributes from a reference object if provided
      if (!is.null(reference)) {
        for (att in intersect(names(attributes(reference)),
                              c("Labels", "Size", "Diag", "Upper"))) {
          attr(d, att) <- attr(reference, att)
        }
      }

    } else {
      # Rectangular: cannot coerce to 'dist', keep as matrix
      d <- dist_matrix
      class(d) <- c("dissimilarity", "matrix")
    }

    return(d)
  }

  #
  # distance <- to_dissimilarity(distance)
  return(distance)

}
