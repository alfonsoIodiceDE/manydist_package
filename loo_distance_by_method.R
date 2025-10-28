loo_distance_by_method = function(df, method = c("baseline","naive","gower","hl","u_dep","u_ind","gudmm","dkss","mod_gower"),dims=2){
 source("R/distance_by_method.r") 
source("R/congruence_coeff.r") 
  
   X=df
  
  full_distance = distance_by_method(X, method = method) |> as.matrix()
  MDS_full <- cmdscale(full_distance, eig=TRUE, k=2)
  
  
  loo_distances <- map(names(X), ~ {
    X_subset <- X %>% select(-all_of(.x))
    distance_by_method(X_subset,method) |> as.matrix() 
  })
  
  
  
  names(loo_distances) <- names(X)
  
 ######### Now I need to add all sorts of computations
  loo_res = tibble(loo_var=names(loo_distances), loo_dist=loo_distances) |> 
    mutate(
      mad_importances = map_dbl(loo_dist, ~ mean(abs(full_distance - .x))),
      cc_importances = map_dbl(loo_dist, ~ {
        MDS_loo <- cmdscale(.x, eig=TRUE, k=dims)
        congruence_coeff(MDS_full$points[,1:dims], MDS_loo$points[,1:dims])
      }),
      ac_importances = sqrt(1 - cc_importances^2),
      mad_normalized = mad_importances / sum(mad_importances)  # like maddis
    ) |> 
    select(-loo_dist)
  
  
  return(loo_res)

}
