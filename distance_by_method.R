distance_by_method <-  function(df, method = c("baseline","naive","gower","hl","u_dep","u_ind","gudmm","dkss","mod_gower")){
  source("R/gudmm_preprocessing.R")
  source("R/gudmm_distance_dependency_mixed_matrix.R")
  source("R/mg_gower_mod_matrix.R")
  source("R/dkss_preprocessing.R")
  
  
  
  
  
  
  no_f_cont <-  df |> dplyr::select(where(is.numeric)) |> ncol()
  
  no_f_ord <- 0 #df |> select(where(is.ordered)) |> ncol()
  
  if (method == "baseline"){
    
  distance_matrix <-mdist(df,distance_cont = "manhattan",commensurable=FALSE, scaling_cont ="std")
  
  }else if (method == "naive"){
    # Naive approach: Euclidean distance on one-hot encoded data
    distance_matrix <- mdist(df, preset = "euclidean_onehot")
    
  }else if (method == "gower") {
    # Standard Gower distance
    distance_matrix <- mdist(df, preset = "gower") * ncol(df)
  } else if (method == "hl") {
    # Hennig-Liao distance with Manhattan for continuous and HL for categorical
    distance_matrix <- mdist(df, 
                             preset = "custom",
                             distance_cont = "euclidean",
                             distance_cat = "HLeucl", 
                             commensurable = FALSE, 
                             scaling_cont = "std")
  }else if (method == "hl_add") {
    # Hennig-Liao distance with Manhattan for continuous and HL for categorical
    distance_matrix <- mdist(df, 
                             preset = "custom",
                             distance_cont = "manhattan",
                             distance_cat = "HL", 
                             commensurable = FALSE, 
                             scaling_cont = "std")
  } else if (method == "u_dep") {
    # Unbiased dependent: PC scores and Total Variation Distance
      distance_matrix <- mdist(df, preset = "unbiased_dependent")
     # distance_matrix <- mdist(df, distance_cat="tot_var_dist", scaling_cont="pc_scores", commensurable = TRUE)
  } else if (method == "u_ind") {
    # Unbiased independent: standardized Manhattan and matching
    distance_matrix <- mdist(df, 
                             preset = "custom",
                             distance_cont = "manhattan",
                             distance_cat = "matching",
                             commensurable = TRUE,
                             scaling_cont = "std")
  }else if (method == "u_dep_manh") {
    # Unbiased independent: standardized Manhattan and matching
    distance_matrix <- mdist(df, 
                             preset = "custom",
                             distance_cont = "manhattan",
                             distance_cat = "tot_var_dist",
                             commensurable = TRUE,
                             scaling_cont = "std")
  }else if (method == "u_dep_eucl") {
    # Unbiased independent: standardized Manhattan and matching
    distance_matrix <- mdist(df, 
                             preset = "custom",
                             distance_cont = "euclidean",
                             distance_cat = "tot_var_dist",
                             commensurable = TRUE,
                             scaling_cont = "std")
  }
  else if (method == "gudmm") {
    X_matrix <- gudmm_preprocessing(df,no_f_cont)
    # print(head(as_tibble(X_matrix)))
    Di = gudmm_distance_dependency_mixed_matrix(X_matrix, no_f_cont, no_f_ord, 'DM5')
    distance_matrix=as.dist(Di)      
  } else if (method == "dkss") {
    
    df = dkss_preprocessing(df, no_f_cont)
    
    # print(df |> head())
    
    dkss_result <- kdml::dkss(df = df,
                              bw = "np",  # Fast bandwidth selection
                              cFUN = "c_gaussian",
                              uFUN = "u_aitken",
                              oFUN = "o_wangvanryzin",
                              stan = TRUE,
                              verbose = FALSE)
    # Di <- tryCatch({
    # print(dkss_result)
    
    distance_matrix = dkss_result$distances
    
    # }, error = function(e) {
    #   warning(paste("DKSS failed:", e$message))
    #   return(NULL)
    # })
    
  } else if (method == "mod_gower") {
    distance_matrix <- mg_gower_mod_matrix(df, use_weights = TRUE) |> as.dist()
  } else {
    stop("Invalid method specified.")
  }
  return(distance_matrix)
}
