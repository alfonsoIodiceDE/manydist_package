## Add a warning message when weight_sys is not user-defined and weights is differnt that all ones
## Add a check when the length of weights vector is different from nvars

ndist <- function(x, commensurable=FALSE, method = "manhattan",
                  scaling = "none" ,ncomp=ncol(x),threshold=NULL,weights=rep(1,ncol(x))){

  .x <- NULL
  by_var_dist_w <- NULL

  rec_x = recipe(~., data = x)

  if(scaling == "std"){

    x = rec_x |> step_scale(all_predictors()) |> prep(training = x) |> bake(new_data=NULL)

  }else if(scaling == "range"){
    x = rec_x |> step_range(all_predictors()) |> prep(training = x) |> bake(new_data=NULL)
  }else if(scaling == "pc_scores" ){
    if(is.null(threshold)){
      x = rec_x |> step_center(all_predictors()) |> step_scale(all_predictors()) |>
        step_pca(all_predictors(),num_comp=ncomp) |> prep(training = x) |> bake(new_data=NULL)
    }else{
      x = rec_x |> step_center(all_predictors()) |> step_scale(all_predictors()) |>
        step_pca(all_predictors(),threshold=threshold) |> prep(training = x) |> bake(new_data=NULL)
    }
  }else if(scaling=="robust"){
    x = x  |>  map(.f=function(x=.x){
      x=(x-median(x))/IQR(x)
      return(x)
    }) |> dplyr::bind_cols()
  }#else{print(is.matrix(x))}


  if(commensurable == TRUE){
    by_var_dist = map(.x=as_tibble(x),.f=function(x=.x){
      b_v_d = daisy(data.frame(x),metric=method) %>% as.matrix()
      b_v_d = b_v_d/mean(b_v_d)
      return(b_v_d)
    })

    distance = Reduce(`+`,by_var_dist)


  }else if (commensurable==FALSE){
    by_var_dist = map(.x=as_tibble(x),.f=function(x=.x){
      b_v_d = daisy(data.frame(x),metric=method, warnBin = FALSE)
      return(b_v_d)
    })

    if (method == "euclidean") {
      by_var_structure = tibble(by_var_dist=by_var_dist,weights=weights) |>
      #  mutate(by_var_dist_w= map(.x=by_var_dist,~.x^2))
        mutate(by_var_dist_w= map2(.x=by_var_dist,.y=weights,~(.x^2)*.y))

      distance = as.matrix(Reduce(`+`,by_var_structure |> pull(by_var_dist_w)))
      distance = sqrt(distance)

    } else {

      by_var_structure = tibble(by_var_dist=by_var_dist,weights=weights) |>
        mutate(by_var_dist_w= map2(.x=by_var_dist,.y=weights,~.x*.y))

      distance = as.matrix(Reduce(`+`,by_var_structure |> pull(by_var_dist_w)))
    }
  }
  # }else if (weight_sys =="none"){
  #  print('none')
  #    distance = as.matrix(daisy(data.frame(x),metric=method))
  #  } else {
  #   print('weightsys should be commensurable, user-defined or none')
  #  }

  return(distance)

}

