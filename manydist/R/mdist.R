mdist <- function(x,new_data=NULL,response=NULL, distance_cont="manhattan", distance_cat="tot_var_dist",
                  commensurable = FALSE,scaling_cont="none",
                  ncomp=ncol(x), threshold = NULL,preset = "custom"){#,prop_nn=0.1, alpha=.5){


  response_quo <- rlang::enquo(response)

  if (!rlang::quo_is_null(response_quo)) {
    # Resolve the selection on the training data
    resp_idx <- tidyselect::eval_select(response_quo, x)

    if (length(resp_idx) != 1L) {
      stop("`response` must select exactly one column.", call. = FALSE)
    }

    resp_name <- names(resp_idx)

    # Extract response vector
    y <- x[[resp_name]]

    # Remove response column from x
    x <- x[, setdiff(colnames(x), resp_name), drop = FALSE]

    # Also drop response column from new_data if present
    if (!is.null(new_data)) {
      if (resp_name %in% colnames(new_data)) {
        new_data <- new_data[, setdiff(colnames(new_data), resp_name), drop = FALSE]
      } else {
        warning("`response` column '", resp_name,
                "' not found in `new_data`; only removed from training `x`.",
                call. = FALSE)
      }
    }
  } else {
    y <- NULL
  }


  ################################################
  ################################################
  validate_x = new_data
  ################################################
  ################################################

  .x = NULL
  a <- NULL
  b <- NULL
  gowdist <- NULL
  cat_data  = x %>% dplyr::select(where(is.factor))
  cont_data = x %>% dplyr::select(where(is.numeric))

  if (ncol(cat_data) == 0)
    cat_data = NULL

  if (ncol(cont_data) == 0)
    cont_data = NULL

  # --- Apply preset defaults for later consistency with continuous only and categorical only ---
  if (preset == "unbiased_dependent") {
    distance_cont <- "manhattan"
    distance_cat  <- "tot_var_dist"
    commensurable <- TRUE
    scaling_cont  <- "pc_scores"
  }

  if (preset == "euclidean_onehot") {
    distance_cont <- "euclidean"
    commensurable <- FALSE
    scaling_cont  <- "std"
  }


  # --- Robust fallbacks for all presets ---
  if (!is.null(cat_data) && distance_cat == "tot_var_dist" && ncol(cat_data) == 1) {
    warning("'tot_var_dist' requires >1 categorical variable. Switching to 'matching'.", call. = FALSE)
    distance_cat <- "matching"
  }

  if (!is.null(cont_data) && scaling_cont == "pc_scores" && ncol(cont_data) == 1) {
    warning("With 1 continuous variable, 'pc_scores' is equivalent to standardization. Switching to scaling_cont='std'.",
            call. = FALSE)
    scaling_cont <- "std"
    # (or scaling_cont <- "none", depending on how ndist defines "std")
  }
  if(!is.null(validate_x)){
    cat_data_val  = validate_x %>% dplyr::select(where(is.factor))
    cont_data_val = validate_x %>% dplyr::select(where(is.numeric))
    if (ncol(cat_data_val) == 0)
      cat_data_val = NULL

    if (ncol(cont_data_val) == 0)
      cont_data_val = NULL
  }


  #### ACTUAL MIXED DATA
  if(!is.null(cont_data) & !is.null(cat_data)){

    if(preset == "gower"){
      if(is.null(validate_x)){
        if (commensurable == FALSE){
          distance_mat <-  as.matrix(daisy(x, metric = "gower"))
        } else {

          gowerlist = x %>% map(~daisy(as_tibble(.x),metric="gower") %>% as.matrix())
          gowerlist = tibble(gowdist = gowerlist) %>% mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)
          distance_mat <- distance_mat
        }
      }else{
        gow_prep <- gower_recipe(data=x) |> prep(training=x)
        x = gow_prep |> bake(new_data = NULL)
        validate_x = gow_prep |> bake(new_data = validate_x)

        if (commensurable == FALSE){
          # distance_mat <-  as.matrix(dist(x, method = "manhattan"))[1:5,1:5]

          gowerlist = map2(
            .x = as_tibble(x),
            .y = as_tibble(validate_x),
            .f = function(x, xnew) {
              x    <- as.numeric(x)
              xnew <- as.numeric(xnew)

              # 1D Manhattan (and 1D Euclidean) distance:
              b_v_d <- abs(outer(xnew, x, "-"))

            }
          )

          # gowerlist = tibble(gowdist = gowerlist) |>
          #   mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))

          distance_mat <- Reduce(`+`, gowerlist)


        }else{
          gowerlist = map2(
            .x = as_tibble(x),
            .y = as_tibble(validate_x),
            .f = function(x, xnew) {
              x    <- as.numeric(x)
              xnew <- as.numeric(xnew)

              # 1D Manhattan (and 1D Euclidean) distance:
              b_v_d <- abs(outer(xnew, x, "-"))

              b_v_d / mean(b_v_d)
            }
          )

          # gowerlist = tibble(gowdist = gowerlist) |>
          #   mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))

          distance_mat <- Reduce(`+`, gowerlist)
        }
      }
    }  else if(preset == "unbiased_dependent"){
      distance_cont = "manhattan"
      distance_cat = "tot_var_dist"
      commensurable = TRUE
      scaling_cont = "pc_scores"


      if(is.null(validate_x)){
        cont_dist_mat = ndist(x = cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,method=distance_cat,commensurable = commensurable)$distance_mat |> as.matrix()

        distance_mat = cat_dist_mat + cont_dist_mat

      }else{

        cont_dist_mat = ndist(x = cont_data, validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,validate_x= cat_data_val,method=distance_cat,commensurable = commensurable)$distance_mat |>  as.matrix()

        distance_mat = cat_dist_mat + cont_dist_mat
      }

      if(interaction){
        i_distance_matrix <- idist(D = cont_dist_mat, cat_data = df_cat, pi_nn = 0.2, decision = "prior_corrected")
        rho_w = 1/ncol(df_cat)
        # distance_matrix = cont_dist_mat + (1-rho_w)*cat_dist_mat + rho_w * i_distance_matrix
        distance_mat = cont_dist_mat + (ncol(df_cat))/(ncol(df_cont) )*(cat_dist_mat + i_distance_matrix)
        # distance_matrix = cont_dist_mat + cat_dist_mat + rho_w*i_distance_matrix
      }

    } else if(preset == "euclidean_onehot"){

      distance_cont = "euclidean"
      commensurable = FALSE
      scaling_cont="std"

      dummy_recipe = recipe(~.,data=cat_data) |> step_dummy(all_nominal(),one_hot = TRUE)

      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data=NULL)

      if(is.null(validate_x)){
        cont_dist_mat = ndist(x=cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = ndist(x=cat_data_dummy, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(x=cont_data,validate_x=cont_data_val,
                              method = distance_cont,
                              commensurable = commensurable,
                              scaling=scaling_cont)  |>  as.matrix()

        cat_data_val_dummy = dummy_recipe |>
          prep(training = cat_data) |>
          bake(new_data=cat_data_val)

        cat_dist_mat = ndist(x=cat_data_dummy,
                             validate_x=cat_data_val_dummy,
                             method = distance_cont,
                             commensurable = commensurable,
                             scaling=scaling_cont)  |>  as.matrix()

        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      }
    }else if(preset=="custom"){
      #    if(weight_cont != "commensurable"){weight_cont = weight_cont}
      #    if(weight_cat != "commensurable"){weight_cat = weight_cat}
      if(is.null(validate_x)){
        cont_dist_mat = ndist(x=cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont, ncomp = ncomp, threshold=threshold)  |>  as.matrix()
        cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat |>  as.matrix()
        distance_mat = cat_dist_mat + cont_dist_mat
        if ((distance_cont == "euclidean") & (distance_cat=="HLeucl")){
          distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        }
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(x=cont_data,validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont, ncomp = ncomp, threshold=threshold)  |>  as.matrix()
        cat_dist_mat = cdist(x=cat_data,validate_x=cat_data_val, response=response,method=distance_cat,commensurable = commensurable)$distance_mat |>  as.matrix()
        distance_mat = cat_dist_mat + cont_dist_mat
        if ((distance_cont == "euclidean") & (distance_cat=="HLeucl"))
          distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      }
      ############################################
      ############ EXTERNAL FUNCTIONS ############
      ############################################

    } else if (preset == "gudmm") {
      no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
      df <- cbind(cont_data, cat_data)
      if (is.null(validate_x)) {
        X_matrix <- gudmm_preprocessing(df, no_f_cont)
        Di <- gudmm_distance_dependency_mixed_matrix(X_matrix, no_f_cont, no_f_ord = 0, DM = "DM5")
        distance_mat <- as.matrix(Di)
      } else {
        stop("train to test distances not implemented for this method")
      }

    } else if (preset == "dkss") {
      no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
      df <- cbind(cont_data, cat_data)
      if (is.null(validate_x)) {
        df <- dkss_preprocessing(df, no_f_cont)
        dkss_result <- kdml::dkss(
          df = df, bw = "np",
          cFUN = "c_gaussian", uFUN = "u_aitken", oFUN = "o_wangvanryzin",
          stan = TRUE, verbose = FALSE
        )
        distance_mat <- dkss_result$distances
      } else {
        stop("train to test distances not implemented for this method")
      }

    } else if (preset == "mod_gower") {
      if (is.null(validate_x)) {
        df <- cbind(cont_data, cat_data)     # no `no_f_cont` needed here
        distance_mat <- mg_gower_mod_matrix(df, use_weights = TRUE) |> as.matrix()
      } else {
        stop("train to test distances not implemented for this method")
      }
    }
    #########


  } else if(is.null(cont_data) & !is.null(cat_data)){### categorical only

    if(preset == "gower"){
      if(is.null(validate_x)){
        if (commensurable == FALSE){
          distance_mat=ncol(cat_data)*daisy(cat_data,metric = "gower") %>% as.matrix()
        }else{
          gowerlist = cat_data %>% map(~daisy(as_tibble(.x),metric="gower") %>% as.matrix())
          gowerlist = tibble(gowdist = gowerlist) %>% mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)
          distance_mat <- distance_mat
        }
      }else{
        gow_prep <- gower_recipe(data=cat_data) |> prep(training=cat_data)
        cat_data = gow_prep |> bake(new_data = NULL)
        cat_data_val = gow_prep |> bake(new_data = cat_data_val)

        if (commensurable == FALSE){

          gowerlist = map2(
            .x = as_tibble(cat_data),
            .y = as_tibble(cat_data_val),
            .f = function(x, xnew) {
              x    <- as.numeric(x)
              xnew <- as.numeric(xnew)

              # 1D Manhattan (and 1D Euclidean) distance:
              b_v_d <- abs(outer(xnew, x, "-"))


            }
          )

          # gowerlist = tibble(gowdist = gowerlist) |>
          #   mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))

          distance_mat <- Reduce(`+`, gowerlist)
        }else{
          # gowerlist = map2(.x=cat_data,.y=cat_data_val,
          #                  ~Rfast::dista(xnew = .y,x = .x,
          #                                type = "manhattan") |> as.matrix()
          # )
          gowerlist = map2(
            .x = as_tibble(cat_data),
            .y = as_tibble(cat_data_val),
            .f = function(x, xnew) {
              x    <- as.numeric(x)
              xnew <- as.numeric(xnew)

              # 1D Manhattan (and 1D Euclidean) distance:
              b_v_d <- abs(outer(xnew, x, "-"))

              b_v_d / mean(b_v_d)
            }
          )

          # gowerlist = tibble(gowdist = gowerlist) |>
          #   mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))

          distance_mat <- Reduce(`+`, gowerlist)

        }
      }

    }else if(preset == "unbiased_dependent"){

      distance_cat = "tot_var_dist"
      commensurable=TRUE
      if(is.null(validate_x)){
        distance_mat = cdist(x = cat_data,method=distance_cat,commensurable = commensurable)$distance_mat |>  as.matrix()

      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        distance_mat = cdist(x = cat_data,validate_x= cat_data_val,method=distance_cat,commensurable = commensurable)$distance_mat |>  as.matrix()
      }

    }else if(preset == "euclidean_onehot"){


      distance_cont <- "euclidean"
      commensurable <- FALSE
      scaling_cont  <- "std"



      dummy_recipe = recipe(~.,data=cat_data) |> step_dummy(all_nominal(),one_hot = TRUE)

      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data=NULL)

      if(is.null(validate_x)){

        distance_mat = ndist(x=cat_data_dummy, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()

      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT


        cat_data_val_dummy = dummy_recipe |>
          prep(training = cat_data) |>
          bake(new_data=cat_data_val)

        distance_mat = ndist(x=cat_data_dummy,
                             validate_x=cat_data_val_dummy,
                             method = distance_cont,
                             commensurable = commensurable,
                             scaling=scaling_cont)  |>  as.matrix()

      }
    }else if(preset=="custom"){

      if(is.null(validate_x)){

        distance_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat |>  as.matrix()


      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT

        distance_mat = cdist(x=cat_data,validate_x=cat_data_val, response=response,method=distance_cat,commensurable = commensurable)$distance_mat |>  as.matrix()


      }
    }else if (preset == "gudmm") {
      no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
      df <- cat_data
      if (is.null(validate_x)) {
        X_matrix <- gudmm_preprocessing(df, no_f_cont)
        Di <- gudmm_distance_dependency_mixed_matrix(X_matrix, no_f_cont, no_f_ord = 0, method = "DM5")
        distance_mat <- as.matrix(Di)
      } else {
        stop("train to test distances not implemented for this method")
      }

    } else if (preset == "dkss") {
      no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
      df <- cat_data
      if (is.null(validate_x)) {
        df <- dkss_preprocessing(df, no_f_cont)
        dkss_result <- kdml::dkss(
          df = df, bw = "np",
          cFUN = "c_gaussian", uFUN = "u_aitken", oFUN = "o_wangvanryzin",
          stan = TRUE, verbose = FALSE
        )
        distance_mat <- dkss_result$distances
      } else {
        stop("train to test distances not implemented for this method")
      }

    } else if (preset == "mod_gower") {
      if (is.null(validate_x)) {
        df <- cat_data
        distance_mat <- mg_gower_mod_matrix(df, use_weights = TRUE) |> as.matrix()
      } else {
        stop("train to test distances not implemented for this method")
      }
    }

    # Continuous only

  }else if(!is.null(cont_data) & is.null(cat_data)){

    if(preset == "gower"){
      if(is.null(validate_x)){
        if (commensurable == FALSE){
          distance_mat <-  ncol(cont_data)*as.matrix(daisy(cont_data, metric = "gower")) # check if this ncol is needed
        } else {

          gowerlist = cont_data %>% map(~daisy(as_tibble(.x),metric="gower") %>% as.matrix())
          gowerlist = tibble(gowdist = gowerlist) %>% mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)

        }
      }else{
        gow_prep <- gower_recipe(data=cont_data) |> prep(training=cont_data)
        cont_data = gow_prep |> bake(new_data = NULL)
        cont_data_val = gow_prep |> bake(new_data = cont_data_val)

        if (commensurable == FALSE){
          # distance_mat <-  as.matrix(dist(x, method = "manhattan"))[1:5,1:5]
          gowerlist = map2(
            .x = as_tibble(cont_data),
            .y = as_tibble(cont_data_val),
            .f = function(x, xnew) {
              x    <- as.numeric(x)
              xnew <- as.numeric(xnew)

              # 1D Manhattan (and 1D Euclidean) distance:
              b_v_d <- abs(outer(xnew, x, "-"))


            }
          )

          # gowerlist = tibble(gowdist = gowerlist) |>
          #   mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))

          distance_mat <- Reduce(`+`, gowerlist)

        }else{
          gowerlist = map2(
            .x = as_tibble(cont_data),
            .y = as_tibble(cont_data_val),
            .f = function(x, xnew) {
              x    <- as.numeric(x)
              xnew <- as.numeric(xnew)

              # 1D Manhattan (and 1D Euclidean) distance:
              b_v_d <- abs(outer(xnew, x, "-"))

              b_v_d / mean(b_v_d)
            }
          )

          # gowerlist = tibble(gowdist = gowerlist) |>
          #   mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))

          distance_mat <- Reduce(`+`, gowerlist)
        }
      }
    }else if(preset == "unbiased_dependent"){
      distance_cont = "manhattan"
      commensurable = TRUE
      scaling_cont = "pc_scores"


      if(is.null(validate_x)){
        distance_mat = ndist(x = cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()

      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        distance_mat = ndist(x = cont_data, validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
      }
    }else if(preset == "euclidean_onehot"){
      distance_cont = "euclidean"
      commensurable = FALSE
      scaling_cont="std"



      if(is.null(validate_x)){
        distance_mat = ndist(x=cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()

      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        distance_mat = ndist(x=cont_data,validate_x=cont_data_val,
                             method = distance_cont,
                             commensurable = commensurable,
                             scaling=scaling_cont,
                             interaction=FALSE)  |>  as.matrix()

      }

    }else if(preset=="custom"){

      if(is.null(validate_x)){
        distance_mat = ndist(x=cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont, ncomp = ncomp, threshold=threshold)  |>  as.matrix()

      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        distance_mat = ndist(x=cont_data,validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont, ncomp = ncomp, threshold=threshold)  |>  as.matrix()
      }

    }else if (preset == "gudmm") {
      no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
      df <- cont_data
      if (is.null(validate_x)) {
        X_matrix <- gudmm_preprocessing(df, no_f_cont)
        Di <- gudmm_distance_dependency_mixed_matrix(X_matrix, no_f_cont, no_f_ord = 0, method = "DM5")
        distance_mat <- as.matrix(Di)
      } else {
        stop("train to test distances not implemented for this method")
      }

    } else if (preset == "dkss") {
      no_f_cont <- if (!is.null(cont_data)) ncol(cont_data) else 0
      df <- cont_data
      if (is.null(validate_x)) {
        df <- dkss_preprocessing(df, no_f_cont)
        dkss_result <- kdml::dkss(
          df = df, bw = "np",
          cFUN = "c_gaussian", uFUN = "u_aitken", oFUN = "o_wangvanryzin",
          stan = TRUE, verbose = FALSE
        )
        distance_mat <- dkss_result$distances
      } else {
        stop("train to test distances not implemented for this method")
      }

    } else if (preset == "mod_gower") {
      if (is.null(validate_x)) {
        df <- cont_data
        distance_mat <- mg_gower_mod_matrix(df, use_weights = TRUE) |> as.matrix()
      } else {
        stop("train to test distances not implemented for this method")
      }
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

  distance_mat <- to_dissimilarity(distance_mat)

  params <- list(
    distance_cont = distance_cont,
    distance_cat  = distance_cat,
    scaling_cont  = scaling_cont,
    commensurable = commensurable,
    ncomp         = ncomp,
    threshold     = threshold,
    rectangular   = (nrow(as.matrix(distance_mat)) != ncol(as.matrix(distance_mat))),
    train_n       = nrow(x),
    test_n        = if (!is.null(validate_x)) nrow(validate_x) else NULL,
    cat_p         = if (!is.null(cat_data))  ncol(cat_data)  else 0,
    cont_p        = if (!is.null(cont_data)) ncol(cont_data) else 0,
    response_col  = if (!rlang::quo_is_null(response_quo)) resp_name else NULL,
    y             = y  # optional: store labels if useful
  )

  return(MDist$new(
    distance = distance_mat,
    preset   = preset,
    data     = x,        # or x if you want to keep it
    params   = params
  ))

}


