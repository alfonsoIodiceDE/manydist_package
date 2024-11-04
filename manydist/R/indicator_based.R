indicator_based<-function(cat_data,commensurable=FALSE, scaling="none", weights=1){
  {
    factor_name = NULL
    delta = NULL
    by_var_dist = NULL
    mean_by_var_dist = NULL
    .x = NULL
    eta = NULL
    scaled_Zs = NULL
    comm_dist = NULL

    #  if (weight_sys != "commensurable") {
    delta_output = cat_delta(x=cat_data,method_cat="matching")
    Z = delta_output$Z
    delta_m = delta_output$matching
    delta_ms <- NULL
    delta_names = delta_output$delta_names
    p <- ncol(cat_data)

    if(scaling=="none"){
      delta_m = 2 * delta_m#sqrt(2) * delta_m
      #  cat_dist_mat = Z %*%  delta_m  %*% t(Z)

      if(!commensurable) {
        cat_dist_mat = Z %*%  delta_m  %*% t(Z)
      } else {
        #### Commensurability
        Z_list = cat_data |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)) |>
            bake(new_data=NULL)
        )

        Q=map_dbl(cat_data,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)

        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta[levels_identifier==.x,
                                    levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )

        cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))

        #       print(cat_dist_mat[1:4,1:4])
      }

    }else if(scaling=="st_dev"){

      # Create the qs_vec
      qs_vec = cat_data %>%
        as_tibble() %>%
        map(function(x) { as.vector(rep(nlevels(x), nlevels(x))) }) %>%
        unlist()

      # Create the qs_diag
      qs_diag = diag(1 / unlist(qs_vec))

      # Dummy encode and center categorical features
      c_Z = cat_data %>%
        as_tibble() %>%
        recipe(~.) %>%
        step_dummy(all_predictors(), one_hot = TRUE) %>%
        step_center(all_predictors()) %>%
        prep(training = cat_data) %>%
        bake(new_data = NULL) %>%
        as.matrix()

      # Calculate the covariance matrix S
      S = (1 / nrow(c_Z)) * t(c_Z) %*% c_Z  # EQ 17

      # Calculate the inverse square root of the diagonal of S
      inv_sq_S_d = diag(1 / sqrt(diag(S)))

      # Dummy encode without centering
      Z = cat_data %>%
        recipe(~.) %>%
        step_dummy(all_nominal(), one_hot = TRUE) %>%
        prep(training = cat_data) %>%
        bake(new_data = NULL) %>%
        as.matrix()

      # Convert Z to numeric
      Z = apply(Z, 2, as.numeric)

      delta_ms <- ncol(cat_data)*inv_sq_S_d %*% qs_diag %*% delta_m %*% inv_sq_S_d
      # Compute the category dissimilarity scaled distance matrix
      #   cat_dist_mat = Z %*% delta_sd %*% t(Z)
      #      print(cat_dist_mat[1:4,1:5])
      if(!commensurable) {
        # Compute the category dissimilarity scaled distance matrix
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
        #    print(cat_dist_mat[1:4,1:5])
      } else {

        #### Commensurability
        Z_list = cat_data |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)) |>
            bake(new_data=NULL)
        )

        Q=map_dbl(cat_data,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)

        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )

        cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
      }

    }else if(scaling=="HL"){

      eta_vec = cat_data |>as_tibble() |>
        map(function(x=.x){
          as.vector(rep(distancefactor(
            cat=nlevels(x),catsizes=table(x)
          ),nlevels(x))
          )
        }
        )

      eta_diag = diag(unlist(eta_vec))
      delta_ms <- 2*delta_m %*% eta_diag
      #     print(delta_m)
      #      print(delta_ms)
      #  cat_dist_mat = Z %*% delta_hl  %*% t(Z)

      if(!commensurable) {
        cat_dist_mat = data.matrix(Z) %*% delta_ms  %*% t(data.matrix(Z))
      } else {
        #### Commensurability
        Z_list = cat_data |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)) |>
            bake(new_data=NULL)
        )

        Q=map_dbl(cat_data,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)

        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )

        cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
      }

    }else if(scaling=="cat_dis"){

      # Create the qs_vec
      qs_vec = cat_data %>%
        as_tibble() %>%
        map(function(x) { as.vector(rep(nlevels(x), nlevels(x))) }) %>%
        unlist()

      # Create the qs_diag
      qs_diag = diag(1 / unlist(qs_vec))

      # Dummy encode and center categorical features
      c_Z = cat_data %>%
        as_tibble() %>%
        recipe(~.) %>%
        step_dummy(all_predictors(), one_hot = TRUE) %>%
        step_center(all_predictors()) %>%
        prep(training = cat_data) %>%
        bake(new_data = NULL) %>%
        as.matrix()

      # Calculate the covariance matrix S
      S = (1 / nrow(c_Z)) * t(c_Z) %*% c_Z  # EQ 17

      # Calculate the inverse square root of the diagonal of S
      inv_sq_S_d = diag(1 / sqrt(diag(S)))

      # Dummy encode without centering
      Z = cat_data %>%
        recipe(~.) %>%
        step_dummy(all_nominal(), one_hot = TRUE) %>%
        prep(training = cat_data) %>%
        bake(new_data = NULL) %>%
        as.matrix()

      # Convert Z to numeric
      Z = apply(Z, 2, as.numeric)

      Zs = Z %*% inv_sq_S_d %*% diag(1/sqrt(unlist(qs_vec))) # EQ 18

      delta_ms = sqrt(1/qs_vec)*(inv_sq_S_d %*% delta_m %*% inv_sq_S_d)
      # Compute the category dissimilarity scaled distance matrix
      # cat_dist_mat = Z %*% inv_sq_S_d %*% delta_m %*% inv_sq_S_d %*% t(Z)
      #      print(cat_dist_mat[1:4,1:4])
      if(!commensurable) {
        # Compute the category dissimilarity scaled distance matrix
        cat_dist_mat = Z %*% delta_ms %*% t(Z)
        #   print(cat_dist_mat[1:4,1:5])
      } else {
        #### Commensurability
        Z_list = cat_data |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)) |>
            bake(new_data=NULL)
        )

        Q=map_dbl(cat_data,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)

        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_ms[levels_identifier==.x,
                                       levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )

        cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))

      }
    }else if(scaling=="HLeucl"){

      delta_ms <- delta_m
      Z <- cat_data %>%
        recipe(~.) %>%
        step_select(all_nominal()) %>%
        step_dummy(all_nominal(), one_hot = TRUE) %>%
        prep() %>%
        bake(new_data = NULL)

      ### PREVIOUS eta_vec 7-Oct-24
      #  eta_vec = cat_data |>as_tibble() |>
      #    map(function(x=.x){
      #      as.vector(rep(distancefactor(
      #        cat=nlevels(x),catsizes = table(x)
      #      ),nlevels(x))
      #      )
      #    }
      #    )

      eta_vec = cat_data |>as_tibble() |>
        map(function(x=.x){
          as.vector(rep(distancefactor(
            cat=nlevels(x),catsizes = table(x)
          ),nlevels(x))
          )
        }
        )
      # Apply Hennig and Liao scaling to categorical data
      Zs <- (as.matrix(Z) %*% diag(unlist(eta_vec)))

      ### PREVIOUS eta_vec 7-Oct-24
      # Zs <- (as.matrix(Z) %*% diag(unlist(eta_vec))/2)

      if(!commensurable) {
        #   if(weight_cat != "commensurable") {
        cat_dist_mat = daisy(Zs,metric = "euclidean") |> as.matrix()
        #    }
      } else {
        ####### NEW CODE

        # Dummy encode categorical variables
        Z_list <- cat_data %>%
          map(~ as_tibble(.x) %>%
                recipe(~.) %>%
                step_dummy(all_nominal(), one_hot = TRUE) %>%
                prep() %>%
                bake(new_data = NULL)
          )

        # Apply Hennig and Liao scaling to each categorical variable
        eta_vec <- cat_data %>%
          as_tibble() %>%
          map(function(x) {
            as.vector(rep(distancefactor(
              cat = nlevels(x),
              catsizes = table(x)
            ), nlevels(x)))
          })

        # Create a distance matrix for each categorical variable
        commensurable_dist_structure <- tibble(factor_name = names(cat_data)) %>%
          mutate(
            Zs = Z_list,
            eta = eta_vec,
         #   scaled_Zs = map2(Zs, eta, ~ as.matrix(.x) %*% diag(.y) / 2), # Hennig-Liao scaling
            scaled_Zs = map2(Zs, eta, ~ as.matrix(.x) %*% diag(.y)), # Hennig-Liao scaling
            by_var_dist = map(scaled_Zs, ~ daisy(.x, metric = "euclidean") %>% as.matrix()), # Compute distance matrix
            mean_by_var_dist = map_dbl(by_var_dist, ~ mean(.x)), # Compute mean distance
            comm_dist = map2(by_var_dist, mean_by_var_dist, ~ .x / .y) # Normalize for commensurability
          )
        # Aggregate the commensurable distance matrices by summing them up
        cat_dist_mat <- Reduce(`+`, commensurable_dist_structure %>% pull(comm_dist))
      }
      ######


    }else if(scaling=="mca"){

      # Create the qs_vec
      qs_vec <- cat_data %>%
        dplyr::select(where(is.factor)) %>%
        as_tibble() %>%
        map(function(x) { as.vector(rep(nlevels(x), nlevels(x))) }) %>%
        unlist()

      # Dummy encode categorical features without centering
      Z <- cat_data %>%
        recipe(~.) %>%
        step_select(all_nominal()) %>%
        step_dummy(all_nominal(), one_hot = TRUE) %>%
        prep(training = cat_data) %>%
        bake(new_data = NULL) %>%
        as.matrix()

      # Convert the matrix to numeric to avoid any issues with factors/characters
      Z <- apply(Z, 2, as.numeric)

      # Calculate the diagonal matrix of category proportions
      D_p <- diag(colSums(Z) / nrow(Z))

      # Calculate the inverse square root of D_p
      inv_sq_D_p <- diag(1 / sqrt(diag(D_p)))
      #      print('INSIDE')
      # Compute the category dissimilarity scaled distance matrix
      #cat_dist_mat <- Z %*% inv_sq_D_p %*% delta_m %*% inv_sq_D_p %*% t(Z)

      if(!commensurable) {
        cat_dist_mat = Z %*% inv_sq_D_p %*% delta_m %*% inv_sq_D_p %*% t(Z)
      } else {
        #### Commensurability
        Z_list = cat_data |> map(
          ~as_tibble(.x) |> recipe(~.)|>
            step_dummy(all_predictors(),one_hot = TRUE) |>
            prep(training = as_tibble(.x)) |>
            bake(new_data=NULL)
        )

        Q=map_dbl(cat_data,nlevels)
        levels_identifier = rep(names(Q), times = as.vector(Q))
        #delta_out = cat_delta(cats,method = method)

        commensurable_dist_structure = tibble(factor_name = names(Q)) |>
          mutate(delta = map(.x=factor_name,
                             ~delta_mca[levels_identifier==.x,
                                        levels_identifier==.x] |> as.matrix()
          ),
          Zs=Z_list,
          by_var_dist = map2(.x=Z_list,.y=delta,
                             ~as.matrix(.x) %*% .y %*% t(as.matrix(.x))),
          mean_by_var_dist = map_dbl(by_var_dist,~mean(.x)),
          comm_dist = map2(.x=by_var_dist,.y=mean_by_var_dist,
                           ~.x /.y)
          )

        cat_dist_mat = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))




        #       print(cat_dist_mat[1:4,1:4])
      }
    }
  }
  #    else {

  #   }

  #}

  out_catdist = list()
  out_catdist$distance_mat = cat_dist_mat
  out_catdist$delta = delta_m
  out_catdist$delta_ms = delta_ms
  out_catdist$delta_names = delta_names
  return(out_catdist)
}
