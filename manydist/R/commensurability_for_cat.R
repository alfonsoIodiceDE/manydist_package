commensurability_for_cat<-function(x, delta) {

#  library(tidymodels)
#  library(tidyverse)
#  library(data.table)
#  library(Matrix)
  .x = NULL
  by_var_dist = NULL
  mean_by_var_dist = NULL
  factor_name = NULL
  by_var_dist = NULL
  comm_dist = NULL
  weight = NULL


  cats = x

  Z_list = cats |> map(
    ~as_tibble(.x) |> recipe(~.)|>
      step_dummy(all_predictors(),one_hot = TRUE) |>
      prep(training = as_tibble(.x)) |>
      bake(new_data=NULL)
  )

  Q=map_dbl(cats,nlevels)
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

  distance = Reduce(`+`,commensurable_dist_structure |> pull(comm_dist))
  return(distance)

}
