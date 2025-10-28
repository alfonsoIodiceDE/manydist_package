mdist_evolved <- function(x,validate_x=NULL,response=NULL, distance_cont="manhattan", distance_cat="tot_var_dist",
                          commensurable = TRUE,scaling_cont="none",
                          ncomp=ncol(x), threshold = NULL,preset = "custom"){#,prop_nn=0.1, alpha=.5){
  
  source("R/gower_recipe.R")
  source("R/ndist.R")
  source("R/cdist.R")
  source("R/gudmm_preprocessing.R")
  source("R/gudmm_distance_dependency_mixed_matrix.R")
  source("R/mg_gower_mod_matrix.R")
  source("R/dkss_preprocessing.R")
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
  
  
  # Check if tot_var_dist is specified but only one categorical variable exists
  if (!is.null(cat_data) && distance_cat == "tot_var_dist" && ncol(cat_data) == 1) {
    warning("'tot_var_dist' requires more than one categorical variable. Switching to 'matching' distance.")
    distance_cat <- "matching"
  }
  
  # Check if pc_scores scaling_cont is specified but only one continuous variable exists
  if (!is.null(cont_data) && scaling_cont == "pc_scores" && ncol(cont_data) == 1) {
    warning("With only one variable, PCA produces a single component identical to standardization. Consider using scaling_cont = \"std\" instead.")
    
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
          
          distance_mat <-  Rfast::dista(xnew = validate_x, 
                                        x = x,
                                        type = "manhattan") |> as.matrix()
          
          
        }else{
          gowerlist = map2(.x=x,.y=validate_x,
                           ~Rfast::dista(xnew = .y,x = .x,
                                         type = "manhattan") |> as.matrix()
          )
          
          gowerlist = tibble(gowdist = gowerlist) |>  
            mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))
          
          distance_mat <- Reduce(`+`, gowerlist$commgow)
          distance_mat <- distance_mat
        }
      }
    }  else if(preset == "unbiased_dependent"){
      distance_cont = "manhattan"
      distance_cat = "tot_var_dist"
      commensurable = TRUE
      scaling_cont = "pc_scores"
      
      
      if(is.null(validate_x)){
        cont_dist_mat = ndist(x = cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(x = cont_data, validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont)  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,validate_x= cat_data_val,method=distance_cat,commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
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
        cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
        if ((distance_cont == "euclidean") & (distance_cat=="HLeucl")){
          distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        }
      }else{### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT    
        cont_dist_mat = ndist(x=cont_data,validate_x=cont_data_val, method = distance_cont,commensurable = commensurable,scaling=scaling_cont, ncomp = ncomp, threshold=threshold)  |>  as.matrix()
        cat_dist_mat = cdist(x=cat_data,validate_x=cat_data_val, response=response,method=distance_cat,commensurable = commensurable)$distance_mat
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
        Di <- gudmm_distance_dependency_mixed_matrix(X_matrix, no_f_cont, no_f_ord = 0, method = "DM5")
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
      
      distance_mat=ncol(cat_data)*daisy(cat_data,metric = "gower") %>% as.matrix()
      
      
    }else if(preset == "unbiased_dependent"){
      
      distance_cat = "tot_var_dist"
      commensurable=TRUE
      cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable=commensurable)$distance_mat
      distance_mat = cat_dist_mat
      
    }else if(preset == "euclidean_onehot"){
      
      dummy_recipe = recipe(~.,data=cat_data) |> step_dummy(all_nominal(),one_hot = TRUE)
      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data=NULL)
      
      distance_mat = ndist(cat_data_dummy, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
      
    }else if(preset=="custom"){
      
      distance_mat = cdist(cat_data, method = distance_cat,commensurable=commensurable)$distance_mat  |>  as.matrix()
      
    }else if(preset %in% c("gudmm","dkss","mod_gower")){
      stop("the selected method applies to  mixed datasets only")
      }
    
    # Continuous only
    
  }else if(!is.null(cont_data) & is.null(cat_data)){
    
    if(preset == "gower"){
      
      distance_mat=ncol(cont_data)*daisy(cont_data,metric= "gower") %>% as.matrix()
      
    }else if(preset == "unbiased_dependent"){
      distance_cont = "manhattan"
      commensurable=TRUE
      scaling_cont="std"
      cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
      distance_mat = cont_dist_mat
    }else if(preset == "euclidean_onehot"){
      distance_cont = "euclidean"
      commensurable=FALSE
      scaling_cont="std"
      distance_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
    
    }else if(preset=="custom"){
      distance_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling_cont)  |>  as.matrix()
      
    }else if(preset %in% c("gudmm","dkss","mod_gower")){
        stop("the selected method applies to  mixed datasets only")
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
  return(distance_mat)
  
}


