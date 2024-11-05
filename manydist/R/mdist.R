mdist <- function(x,response=NULL, distance_cont="manhattan", distance_cat="tot_var_dist",
                  commensurable = FALSE,scaling="none",
                  ncomp=ncol(x), threshold = NULL,preset = "custom"){#,prop_nn=0.1, alpha=.5){

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

  #### ACTUAL MIXED DATA
  if(!is.null(cont_data) & !is.null(cat_data)){
    if(preset == "gower"){
      if (commensurable == FALSE)
      {
        distance_mat <-  as.matrix(daisy(x, metric = "gower"))
      } else {
        gowerlist = x %>% map(~daisy(as_tibble(.x),metric="gower") %>% as.matrix())
        gowerlist = tibble(gowdist = gowerlist) %>% mutate(commgow = map(.x=gowdist, ~.x / mean(.x)))
        distance_mat <- Reduce(`+`, gowerlist$commgow)
        distance_mat <- distance_mat
      }

    }else if(preset == "catdissim"){
      distance_cont = "manhattan"
      distance_cat = "matching"
      commensurable = TRUE
      # scaling="none"
      scaling="std"
      cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling)  |>  as.matrix()
      cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
      distance_mat = cat_dist_mat + cont_dist_mat

    }else if(preset == "unbiased_dependent"){
      distance_cont = "manhattan"
      distance_cat = "tot_var_dist"
      commensurable = TRUE
      # cont_scaling = "pc_scores"
      # cont_scaling="none"
     # cat_scaling="none"
      scaling="std"
      cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling)  |>  as.matrix()
      cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
      distance_mat = cat_dist_mat + cont_dist_mat

    } else if(preset == "euclidean_onehot"){

      distance_cont = "euclidean"
      commensurable = FALSE
      scaling="std"

      dummy_recipe = recipe(~.,data=cat_data) |> step_dummy(all_nominal(),one_hot = TRUE)

      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data=NULL)
      #   cont_dist_mat = ndist(cont_data, method = distance_cont,weight_sys = weight_cont,scaling=cont_scaling)  |>  as.matrix()
      #    cat_dist_mat = ndist(cat_data_dummy, method = distance_cont,weight_sys = weight_cont,scaling=cont_scaling)  |>  as.matrix()
      #    distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))

      cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling)  |>  as.matrix()
      cat_dist_mat = ndist(cat_data_dummy, method = distance_cont,commensurable = commensurable,scaling=scaling)  |>  as.matrix()
      distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
#      print(sum(distance_mat2-distance_mat))
      #   print(sum(joint_distance_mat - distance_mat))
      #  print(all_equal(joint_distance_mat,distance_mat))

    }else if(preset=="entropy_based"){
      # if(weight_cont != "commensurable"){weight_cont = weight_cont}
      #  if(weight_cat != "commensurable"){weight_cat = weight_cat}

      n_cont=ncol(cont_data)
      x=cbind(cont_data,cat_data)
     # distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()

      #   }else if(preset=="indicator_based"){
      #      if(distance_cont!="manhattan"){distance_cont = distance_cont}
      #      if(distance_cat!="tot_var_dist"){distance_cat = distance_cat}
      #  if(weight_cont != "commensurable"){weight_cont = weight_cont}
      #  if(weight_cat != "commensurable"){weight_cat = weight_cat}
      #     if(cont_scaling!="none"){cont_scaling=cont_scaling}

      #      cat_dist_mat <- indicator_based(x,commensurable = commensurable, scaling=cat_scaling, weights=1)$distance_mat
      #     cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable = commensurable,scaling=cont_scaling)  |>  as.matrix()
      # #  print(cont_dist_mat[1:5,1:5])
      #    if ((distance_cont == "euclidean") | (distance_cat == "euclidean"))
      #      distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      #    else
      #      distance_mat = cat_dist_mat + cont_dist_mat


    }else if (preset=="custom"){
      #    if(weight_cont != "commensurable"){weight_cont = weight_cont}
      #    if(weight_cat != "commensurable"){weight_cat = weight_cat}
      cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable = commensurable,scaling=scaling, ncomp = ncomp, threshold=threshold)  |>  as.matrix()
      cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable = commensurable)$distance_mat
      distance_mat = cat_dist_mat + cont_dist_mat
      if ((distance_cont == "euclidean") & (distance_cat=="HLeucl"))
        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
    }

    #### CATEGORICAL ONLY
  } else if(is.null(cont_data) & !is.null(cat_data)){

    if(preset == "gower"){

      distance_mat=ncol(cat_data)*daisy(cat_data,metric = "gower") %>% as.matrix()


    }else if(preset == "catdissim"){

      distance_cat = "matching"
      commensurable = TRUE
      cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable=commensurable)$distance_mat
      distance_mat = cat_dist_mat

    }else if(preset == "unbiased_dependent"){

      distance_cat = "tot_var_dist"
      weight_cat = "commensurable"
      cat_dist_mat = cdist(x=cat_data,method=distance_cat,commensurable=commensurable)$distance_mat
      distance_mat = cat_dist_mat

    }else if(preset == "euclidean_onehot"){

      dummy_recipe = recipe(~.,data=cat_data) |> step_dummy(all_nominal())
      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data=NULL)

      cat_dist_mat = ndist(cat_data_dummy, method = distance_cont,commensurable=commensurable,scaling=scaling)  |>  as.matrix()
      distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))

    }else if(preset=="entropy_based"){

      n_cont=0
      x = cat_data
     # distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()

    }
    else if(preset=="custom"){

      distance_mat = cdist(cat_data, method = distance_cat,commensurable=commensurable)$distance_mat  |>  as.matrix()

    }
    # Continuous only
  }else if(!is.null(cont_data) & is.null(cat_data)){

    if(preset == "gower"){

      distance_mat=ncol(cont_data)*daisy(cont_data,metric= "gower") %>% as.matrix()

    }else if(preset == "catdissim"){
      distance_cont = "manhattan"
      commensurable=TRUE
      # scaling="none"
      scaling="std"
      cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling)  |>  as.matrix()
      distance_mat = cont_dist_mat

    }else if(preset == "unbiased_dependent"){
      distance_cont = "manhattan"
      commensurable=TRUE
      scaling="std"
      cont_dist_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling)  |>  as.matrix()
      distance_mat = cont_dist_mat
    }else if(preset == "euclidean_onehot"){
      distance_cont = "euclidean"
      commensurable=FALSE
      scaling="std"
      distance_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling)  |>  as.matrix()
    }else if(preset=="entropy_based"){

      n_cont=ncol(cont_data)
      x=cont_data
      #    distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()
    }else if(preset=="custom"){
      distance_mat = ndist(cont_data, method = distance_cont,commensurable=commensurable,scaling=scaling)  |>  as.matrix()

    }

  }

  return(distance_mat)
}

