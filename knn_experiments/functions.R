# to produce a distance matrix which consider the association of x's with y (LDA type) :: similar to supervised_full for categorical data
dann_metric <- function(data, groups) {
  G <-  tibble(groups = groups) %>% recipe( ~ groups) %>%
    step_dummy(groups, one_hot = T) %>%
    prep() %>% bake(new_data = NULL) %>% as.matrix()
  row_sums <- rowsum(data, group = groups)
  group_sizes <- as.vector(table(groups))
  M <- as.matrix(row_sums / group_sizes)
  xbar <- colMeans(data)
  one_x <- matrix(1, nrow(data), 1) %*% xbar
  GM <- G %*% M
  g <- M %>% nrow()
  x_mat <- data |> as.matrix() - GM
  W <- (t(x_mat) %*% x_mat) / (nrow(data) - g)
  B <- (t(GM - one_x) %*% (GM - one_x)) / (g - 1)
  eigen_W <- eigen(W)
  W_inv_sqrt <- eigen_W$vectors %*% diag(1 / sqrt(eigen_W$values)) %*% t(eigen_W$vectors)
  B_s <- W_inv_sqrt %*% B %*% W_inv_sqrt
  sigma <- W_inv_sqrt %*% B_s %*% W_inv_sqrt
  return(sigma)
}

# to produce a distance matrix which consider the association of each x (separately) with y (LDA type) :: similar to supervised for categorical data
dann_metric_sv <- function(variable, groups) {
  G <-  tibble(groups = groups) %>% recipe( ~ groups) %>%
    step_dummy(groups, one_hot = T) %>%
    prep() %>% bake(new_data = NULL) %>% as.matrix()
  M <- tapply(variable, groups, mean)
  xbar <- mean(variable)
  one_x <- rep(xbar, length(variable))
  GM <- G %*% M
  g <- M %>% nrow()
  x_mat <- variable |> as.matrix() - GM
  W <- (t(x_mat) %*% x_mat) / (length(variable) - g)
  B <- (t(GM - one_x) %*% (GM - one_x)) / (g - 1)
  sigma <- B / W
  return(sigma)
}

# to compute the Mahalanobis train to test distance matrix with a custom sigma
maha <- function(train,
                 test,
                 center = TRUE,
                 sig = sig) {
  train <- scale(as.matrix(train), center = center, scale = TRUE)
  test <- scale(as.matrix(test), center = center, scale = TRUE)
  dist <- matrix(0, nrow = nrow(train), ncol = nrow(test))
  for (i in 1:nrow(train)) {
    for (j in 1:nrow(test)) {
      diff <- train[i, ] - test[j, ]
      dist[i, j] <- sqrt(t(diff) %*% sig %*% diff)
    }
  }
  return(dist)
}

# to compute the Mahanattan train to test distance matrix with a custom sigma
mana <- function(train,
                 test,
                 center = TRUE,
                 sig = sig) {
  train <- scale(as.matrix(train), center = center, scale = TRUE)
  test <- scale(as.matrix(test), center = center, scale = TRUE)
  dist <- matrix(0, nrow = nrow(train), ncol = nrow(test))
  for (i in 1:nrow(train)) {
    for (j in 1:nrow(test)) {
      diff <- as.matrix(abs(train[i, ] - test[j, ]))
      dist[i, j] <- sqrt(diag(sig)) %*% diff
    }
  }
  return(dist)
}

# to compute knn predictions (with distance matrix as input)
knn_dist <- function(dist_matrix, train_labels, k) {
  n_test <- ncol(dist_matrix)
  predictions <- vector("character", n_test)
  for (i in 1:n_test) {
    dist_col <- dist_matrix[, i]
    nn_indices <- order(dist_col)[1:k]
    nn_labels <- train_labels[nn_indices]
    predictions[i] <- names(sort(table(nn_labels), decreasing = TRUE))[1]
  }
  return(predictions)
}

#generate mix data with a given association and correlation level
genmixdata <- function(n, 
                       num_continuous,         
                       target_corr = 0.8,         
                       y,
                       num_noise = 0,
                       noise_strength = 0.95,
                       num_categorical = 0,
                       levels_per_cat = 3,
                       cat_assoc_strength = 0.8,
                       num_categorical_noise = 0) {
  
  y <- as.factor(y)
  levels_y <- levels(y)
  num_y_levels <- length(levels_y)
  
  df <- data.frame(row_id = seq_len(n))
  
  # Continuous vars associated with y
  for (i in 1:num_continuous) {
    tc <- target_corr
    group_means <- seq(-1, 1, length.out = num_y_levels) * tc * 2
    x <- numeric(n)
    for (j in 1:num_y_levels) {
      idx <- which(y == levels_y[j])
      x[idx] <- rnorm(length(idx), mean = group_means[j], sd = sqrt(1 - tc^2))
    }
    df[[paste0("ContVar_", i)]] <- x
  }
  
  # Continuous noise vars: orthogonal to y, with tunable variance (i.e., influence)
  if (num_noise > 0) {
    for (k in 1:num_noise) {
      z_raw <- rnorm(n)
      z_orth <- lm(z_raw ~ y)$residuals
      z_scaled <- scale(z_orth) * noise_strength  # controls influence in distance
      df[[paste0("NoiseVar_", k)]] <- z_scaled
    }
  }
  
  # Categorical vars associated with y
  if (num_categorical > 0) {
    for (c in 1:num_categorical) {
      cat_var <- character(n)
      for (j in 1:num_y_levels) {
        idx <- which(y == levels_y[j])
        base_probs <- rep((1 - cat_assoc_strength) / (levels_per_cat - 1), levels_per_cat)
        main_level <- sample(1:levels_per_cat, 1)
        base_probs[main_level] <- cat_assoc_strength
        cat_var[idx] <- sample(
          paste0("Cat", c, "_L", 1:levels_per_cat),
          size = length(idx),
          replace = TRUE,
          prob = base_probs
        )
      }
      df[[paste0("CatVar_", c)]] <- as.factor(cat_var)
    }
  }
  
  # Categorical noise vars: unrelated to y
  if (num_categorical_noise > 0) {
    for (c in 1:num_categorical_noise) {
      cat_noise <- sample(
        paste0("CatNoise", c, "_L", 1:levels_per_cat),
        size = n,
        replace = TRUE
      )
      df[[paste0("CatNoiseVar_", c)]] <- as.factor(cat_noise)
    }
  }
  
  df$row_id <- NULL
  return(df)
}

# to compute the accuracy
compute_accuracy <- function(predicted, actual, show_confusion = FALSE) {
  predicted <- factor(predicted)
  actual <- factor(actual, levels = levels(predicted))
  acc <- mean(predicted == actual)
  if (show_confusion) {
    print(table(Predicted = predicted, Actual = actual))
  }
  return(acc)
}

#alfonso's
cran_mdist <- function(x,
                       validate_x = NULL,
                       response = NULL,
                       distance_cont = "manhattan",
                       distance_cat = "tot_var_dist",
                       commensurable = FALSE,
                       scaling = "none",
                       ncomp = ncol(x),
                       threshold = NULL,
                       preset = "custom") {
  #,prop_nn=0.1, alpha=.5){
  
  # source("gower_recipe.R")
  # source("ndist.R")
  # source("cdist.R")
  
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
  
  
  if (!is.null(validate_x)) {
    cat_data_val  = validate_x %>% dplyr::select(where(is.factor))
    cont_data_val = validate_x %>% dplyr::select(where(is.numeric))
    if (ncol(cat_data_val) == 0)
      cat_data_val = NULL
    
    if (ncol(cont_data_val) == 0)
      cont_data_val = NULL
  }
  
  
  #### ACTUAL MIXED DATA
  if (!is.null(cont_data) & !is.null(cat_data)) {
    if (preset == "gower") {
      if (is.null(validate_x)) {
        if (commensurable == FALSE) {
          distance_mat <-  as.matrix(daisy(x, metric = "gower"))
        } else {
          gowerlist = x %>% map( ~ daisy(as_tibble(.x), metric = "gower") %>% as.matrix())
          gowerlist = tibble(gowdist = gowerlist) %>% mutate(commgow = map(.x =
                                                                             gowdist, ~ .x / mean(.x)))
          distance_mat <- Reduce(`+`, gowerlist$commgow)
          distance_mat <- distance_mat
        }
      } else{
        ### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT : CHECK IF DISTA with gower works
        gow_prep <- cran_gower_recipe(data = x) |> prep(training = x)
        x = gow_prep |> bake(new_data = NULL)
        validate_x = gow_prep |> bake(new_data = validate_x)
        
        if (commensurable == FALSE) {
          # distance_mat <-  as.matrix(dist(x, method = "manhattan"))[1:5,1:5]
          
          distance_mat <-  Rfast::dista(xnew = validate_x,
                                        x = x,
                                        type = "manhattan") |> as.matrix()
          
          
        } else{
          gowerlist = map2(
            .x = x,
            .y = validate_x,
            ~ Rfast::dista(
              xnew = .y,
              x = .x,
              type = "manhattan"
            ) |> as.matrix()
          )
          
          gowerlist = tibble(gowdist = gowerlist) |>
            mutate(commgow = map(.x = gowdist, ~ .x / mean(.x)))
          
          distance_mat <- Reduce(`+`, gowerlist$commgow)
          distance_mat <- distance_mat
        }
      }
      
    } else if (preset == "catdissim") {
      distance_cont = "manhattan"
      distance_cat = "matching"
      commensurable = TRUE
      # scaling="none"
      scaling = "std"
      
      if (is.null(validate_x)) {
        cont_dist_mat = ndist(
          cont_data,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,
                             method = distance_cat,
                             commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      } else{
        ### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(
          x = cont_data,
          validate_x = cont_data_val,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        cat_dist_mat = cdist(
          x = cat_data,
          validate_x = cat_data_val,
          method = distance_cat,
          commensurable = commensurable
        )$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      }
    } else if (preset == "unbiased_dependent") {
      distance_cont = "manhattan"
      distance_cat = "tot_var_dist"
      commensurable = TRUE
      cont_scaling = "pc_scores"
      # cont_scaling="none"
      # cat_scaling="none"
      # scaling="std"
      if (is.null(validate_x)) {
        cont_dist_mat = ndist(
          x = cont_data,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,
                             method = distance_cat,
                             commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      } else{
        ### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(
          x = cont_data,
          validate_x = cont_data_val,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        cat_dist_mat = cdist(
          x = cat_data,
          validate_x = cat_data_val,
          method = distance_cat,
          commensurable = commensurable
        )$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
      }
    } else if (preset == "euclidean_onehot") {
      distance_cont = "euclidean"
      commensurable = FALSE
      scaling = "std"
      
      dummy_recipe = recipe( ~ ., data = cat_data) |> step_dummy(all_nominal(), one_hot = TRUE)
      
      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data = NULL)
      
      if (is.null(validate_x)) {
        cont_dist_mat = ndist(
          x = cont_data,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        cat_dist_mat = ndist(
          x = cat_data_dummy,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      } else{
        ### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(
          x = cont_data,
          validate_x = cont_data_val,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        
        cat_data_val_dummy = dummy_recipe |>
          prep(training = cat_data) |>
          bake(new_data = cat_data_val)
        
        cat_dist_mat = ndist(
          x = cat_data_dummy,
          validate_x = cat_data_val_dummy,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling
        )  |>  as.matrix()
        
        distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        
      }
    } else if (preset == "entropy_based") {
      # if(weight_cont != "commensurable"){weight_cont = weight_cont}
      #  if(weight_cat != "commensurable"){weight_cat = weight_cat}
      #
      # n_cont=ncol(cont_data)
      # x=cbind(cont_data,cat_data)
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
      
      
    } else if (preset == "custom") {
      #    if(weight_cont != "commensurable"){weight_cont = weight_cont}
      #    if(weight_cat != "commensurable"){weight_cat = weight_cat}
      if (is.null(validate_x)) {
        cont_dist_mat = ndist(
          x = cont_data,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling,
          ncomp = ncomp,
          threshold = threshold
        )  |>  as.matrix()
        cat_dist_mat = cdist(x = cat_data,
                             method = distance_cat,
                             commensurable = commensurable)$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
        if ((distance_cont == "euclidean") &
            (distance_cat == "HLeucl")) {
          distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
        }
      } else{
        ### MODIFY TO TAKE THE ASSESSMENT INTO ACCOUNT
        cont_dist_mat = ndist(
          x = cont_data,
          validate_x = cont_data_val,
          method = distance_cont,
          commensurable = commensurable,
          scaling = scaling,
          ncomp = ncomp,
          threshold = threshold
        )  |>  as.matrix()
        cat_dist_mat = cdist(
          x = cat_data,
          validate_x = cat_data_val,
          response = response,
          method = distance_cat,
          commensurable = commensurable
        )$distance_mat
        distance_mat = cat_dist_mat + cont_dist_mat
        if ((distance_cont == "euclidean") &
            (distance_cat == "HLeucl"))
          distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      }
      
      ####
    }
    
  } else if (is.null(cont_data) & !is.null(cat_data)) {
    if (preset == "gower") {
      distance_mat = ncol(cat_data) * daisy(cat_data, metric = "gower") %>% as.matrix()
      
      
    } else if (preset == "catdissim") {
      distance_cat = "matching"
      commensurable = TRUE
      cat_dist_mat = cdist(x = cat_data,
                           method = distance_cat,
                           commensurable = commensurable)$distance_mat
      distance_mat = cat_dist_mat
      
    } else if (preset == "unbiased_dependent") {
      distance_cat = "tot_var_dist"
      weight_cat = "commensurable"
      cat_dist_mat = cdist(x = cat_data,
                           method = distance_cat,
                           commensurable = commensurable)$distance_mat
      distance_mat = cat_dist_mat
      
    } else if (preset == "euclidean_onehot") {
      dummy_recipe = recipe( ~ ., data = cat_data) |> step_dummy(all_nominal())
      cat_data_dummy = dummy_recipe |>
        prep(training = cat_data) |>
        bake(new_data = NULL)
      
      cat_dist_mat = ndist(
        cat_data_dummy,
        method = distance_cont,
        commensurable = commensurable,
        scaling = scaling
      )  |>  as.matrix()
      distance_mat = sqrt((cat_dist_mat^2) + (cont_dist_mat^2))
      
    } else if (preset == "entropy_based") {
      n_cont = 0
      x = cat_data
      # distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()
      
    }
    else if (preset == "custom") {
      distance_mat = cdist(cat_data, method = distance_cat, commensurable = commensurable)$distance_mat  |>  as.matrix()
      
    }
    # Continuous only
    
    
  } else if (!is.null(cont_data) & is.null(cat_data)) {
    if (preset == "gower") {
      distance_mat = ncol(cont_data) * daisy(cont_data, metric = "gower") %>% as.matrix()
      
    } else if (preset == "catdissim") {
      distance_cont = "manhattan"
      commensurable = TRUE
      # scaling="none"
      scaling = "std"
      cont_dist_mat = ndist(
        cont_data,
        method = distance_cont,
        commensurable = commensurable,
        scaling = scaling
      )  |>  as.matrix()
      distance_mat = cont_dist_mat
      
    } else if (preset == "unbiased_dependent") {
      distance_cont = "manhattan"
      commensurable = TRUE
      scaling = "std"
      cont_dist_mat = ndist(
        cont_data,
        method = distance_cont,
        commensurable = commensurable,
        scaling = scaling
      )  |>  as.matrix()
      distance_mat = cont_dist_mat
    } else if (preset == "euclidean_onehot") {
      distance_cont = "euclidean"
      commensurable = FALSE
      scaling = "std"
      distance_mat = ndist(
        cont_data,
        method = distance_cont,
        commensurable = commensurable,
        scaling = scaling
      )  |>  as.matrix()
    } else if (preset == "entropy_based") {
      n_cont = ncol(cont_data)
      x = cont_data
      #    distance_mat=GUDMM(X=x, no_f_cont=n_cont, no_f_ord=0) |> as.matrix()
    } else if (preset == "custom") {
      distance_mat = ndist(
        cont_data,
        method = distance_cont,
        commensurable = commensurable,
        scaling = scaling
      )  |>  as.matrix()
      
    }
    
  }
  
  return(distance_mat)
}

cran_gower_recipe <- function(data) {
  n_cat <- data |> dplyr::select(where(is.factor)) |> ncol()
  n_con <- data |> dplyr::select(where(is.numeric)) |> ncol()
  
  if (n_cat > 0) {
    gow_recipe <- recipe( ~ ., data = data) |>
      step_range(all_numeric_predictors(), min = 0, max = 1) |>
      step_dummy(all_nominal_predictors(),
                 one_hot = TRUE,
                 id = "dummy") |>
      step_mutate(dplyr::across(matches("^f.*_"), ~ .x / 2))
  } else {
    gow_recipe <- recipe( ~ ., data = data) |>
      step_range(all_numeric_predictors(), min = 0, max = 1)
  }
  
  return(gow_recipe)
}


