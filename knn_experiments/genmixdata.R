genmixdata <- function(n,
                       num_continuous,
                       target_corr = 0.8,          # scalar OR length num_continuous
                       y,
                       num_noise = 0,
                       noise_strength = 0.95,
                       num_categorical = 0,
                       levels_per_cat = 3,
                       cat_assoc_strength = 0.8,   # scalar OR length num_categorical
                       cat_redundancy = FALSE,
                       cat_redundancy_noise_prob = 0.15,
                       num_categorical_noise = 0) {

  y <- as.factor(y)
  levels_y <- levels(y)
  num_y_levels <- length(levels_y)

  df <- data.frame(row_id = seq_len(n))

  ## -----------------------------
  ## Handle per-variable target_corr
  ## -----------------------------
  if (length(target_corr) == 1L) {
    target_corr_vec <- rep(target_corr, num_continuous)
  } else {
    if (length(target_corr) != num_continuous) {
      stop("`target_corr` must be either length 1 or length `num_continuous`.")
    }
    target_corr_vec <- target_corr
  }

  ## -----------------------------
  ## Continuous vars associated with y
  ## (heterogeneous variance / shape)
  ## -----------------------------
  for (i in seq_len(num_continuous)) {
    tc <- target_corr_vec[i]
    tc2 <- min(tc^2, 0.999)  # small safety

    group_means <- seq(-1, 1, length.out = num_y_levels) * tc * 2
    x_latent <- numeric(n)

    # Latent variable linearly associated with y
    for (j in seq_len(num_y_levels)) {
      idx <- which(y == levels_y[j])
      x_latent[idx] <- rnorm(length(idx), mean = group_means[j], sd = sqrt(1 - tc2))
    }

    # shape / skewness pattern
    shape_type <- c("normal",
                    "skewed_right",
                    "skewed_left",
                    "heavy_tail")[(i - 1L) %% 4L + 1L]

    # Random overall scale (variance heterogeneity)
    scale_i <- runif(1, 0.5, 3)

    # transform to introduce different shapes
    z <- as.numeric(scale(x_latent))  # standardized latent

    if (shape_type == "normal") {
      x_final <- z * scale_i
    } else if (shape_type == "skewed_right") {
      tmp <- exp(z)
      tmp <- as.numeric(scale(tmp))
      x_final <- tmp * scale_i
    } else if (shape_type == "skewed_left") {
      tmp <- exp(-z)
      tmp <- as.numeric(scale(tmp))
      x_final <- tmp * scale_i
    } else if (shape_type == "heavy_tail") {
      t_noise <- rt(n, df = 3)
      tmp <- z + t_noise
      tmp <- as.numeric(scale(tmp))
      x_final <- tmp * scale_i
    }

    df[[paste0("ContVar_", i)]] <- x_final
  }

  ## -----------------------------------------------
  ## Continuous noise vars: orthogonal to y,
  ## with heterogeneous variance / shape
  ## -----------------------------------------------
  if (num_noise > 0) {
    for (k in seq_len(num_noise)) {
      z_raw  <- rnorm(n)
      z_orth <- lm(z_raw ~ y)$residuals    # orthogonal to y
      z      <- as.numeric(scale(z_orth))

      shape_type <- c("normal",
                      "skewed_right",
                      "skewed_left",
                      "heavy_tail")[(k - 1L) %% 4L + 1L]

      # base scale controlled by noise_strength, plus random variability
      scale_k <- runif(1, 0.5, 1.5) * noise_strength

      if (shape_type == "normal") {
        z_final <- z * scale_k
      } else if (shape_type == "skewed_right") {
        tmp <- exp(z)
        tmp <- as.numeric(scale(tmp))
        z_final <- tmp * scale_k
      } else if (shape_type == "skewed_left") {
        tmp <- exp(-z)
        tmp <- as.numeric(scale(tmp))
        z_final <- tmp * scale_k
      } else if (shape_type == "heavy_tail") {
        t_noise <- rt(n, df = 3)
        tmp <- z + t_noise
        tmp <- as.numeric(scale(tmp))
        z_final <- tmp * scale_k
      }

      df[[paste0("NoiseVar_", k)]] <- z_final
    }
  }

  ## -----------------------------
  ## Handle per-variable cat_assoc_strength
  ## -----------------------------
  if (num_categorical > 0) {
    if (length(cat_assoc_strength) == 1L) {
      cat_assoc_vec <- rep(cat_assoc_strength, num_categorical)
    } else {
      if (length(cat_assoc_strength) != num_categorical) {
        stop("`cat_assoc_strength` must be length 1 or length `num_categorical`.")
      }
      cat_assoc_vec <- cat_assoc_strength
    }
  } else {
    cat_assoc_vec <- numeric(0L)
  }

  ## -----------------------------
  ## Categorical vars associated with y
  ##  - if cat_redundancy = FALSE:
  ##      each CatVar_c has its own strength cat_assoc_vec[c]
  ##  - if cat_redundancy = TRUE:
  ##      CatVar_1 is generated from y,
  ##      CatVar_2..CatVar_num_categorical are noisy copies
  ## -----------------------------
  if (num_categorical > 0) {

    if (!cat_redundancy) {
      ## independent cats with possibly different strengths
      for (c in seq_len(num_categorical)) {
        strength_c <- cat_assoc_vec[c]
        cat_var <- character(n)

        for (j in seq_len(num_y_levels)) {
          idx <- which(y == levels_y[j])
          base_probs <- rep((1 - strength_c) / (levels_per_cat - 1L), levels_per_cat)
          main_level <- sample(seq_len(levels_per_cat), 1L)
          base_probs[main_level] <- strength_c

          cat_var[idx] <- sample(
            paste0("Cat", c, "_L", seq_len(levels_per_cat)),
            size = length(idx),
            replace = TRUE,
            prob = base_probs
          )
        }
        df[[paste0("CatVar_", c)]] <- factor(cat_var)
      }

    } else {
      ## one "base" cat associated with y, others are noisy copies (redundant block)
      strength_1 <- cat_assoc_vec[1L]
      base_cat <- character(n)

      for (j in seq_len(num_y_levels)) {
        idx <- which(y == levels_y[j])
        base_probs <- rep((1 - strength_1) / (levels_per_cat - 1L), levels_per_cat)
        main_level <- sample(seq_len(levels_per_cat), 1L)
        base_probs[main_level] <- strength_1

        base_cat[idx] <- sample(
          paste0("Cat1_L", seq_len(levels_per_cat)),
          size = length(idx),
          replace = TRUE,
          prob = base_probs
        )
      }
      base_cat <- factor(base_cat)
      df[["CatVar_1"]] <- base_cat

      if (num_categorical > 1) {
        for (c in 2:num_categorical) {
          tmp <- base_cat
          flip_idx <- sample(seq_len(n), size = round(cat_redundancy_noise_prob * n))
          levs <- levels(tmp)

          for (i in flip_idx) {
            tmp[i] <- sample(levs[levs != tmp[i]], 1L)
          }
          df[[paste0("CatVar_", c)]] <- droplevels(tmp)
        }
      }
    }
  }

  ## -----------------------------
  ## Categorical noise vars: unrelated to y
  ## -----------------------------
  if (num_categorical_noise > 0) {
    for (c in seq_len(num_categorical_noise)) {
      cat_noise <- sample(
        paste0("CatNoise", c, "_L", seq_len(levels_per_cat)),
        size = n,
        replace = TRUE
      )
      df[[paste0("CatNoiseVar_", c)]] <- factor(cat_noise)
    }
  }

  df$row_id <- NULL
  return(df)
}

# # to compute the accuracy
# compute_accuracy <- function(predicted, actual, show_confusion = FALSE) {
#   predicted <- factor(predicted)
#   actual <- factor(actual, levels = levels(predicted))
#   acc <- mean(predicted == actual)
#   if (show_confusion) {
#     print(table(Predicted = predicted, Actual = actual))
#   }
#   return(acc)
# }
#
# # to compute knn predictions (with distance matrix as input)
# knn_dist <- function(dist_matrix, train_labels, k) {
#   n_test <- ncol(dist_matrix)
#   predictions <- vector("character", n_test)
#   for (i in 1:n_test) {
#     dist_col <- dist_matrix[, i]
#     nn_indices <- order(dist_col)[1:k]
#     nn_labels <- train_labels[nn_indices]
#     predictions[i] <- names(sort(table(nn_labels), decreasing = TRUE))[1]
#   }
#   return(predictions)
# }
#
# make_sim_data <- function(n,
#                           num_signal_cont    = 2,
#                           num_noise_cont     = 8,
#                           num_redundant_cat  = 4,   # + 1 base = 5 cat signal vars
#                           num_cat_noise      = 20,
#                           levels_per_cat     = 3,
#                           cont_signal_sd     = 0.8,
#                           cont_noise_sd      = 1,
#                           cat_assoc_strength = 0.9,
#                           flip_prob_redundant = 0.1) {
#   # y: balanced 3-class factor
#   y <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
#   k <- nlevels(y)
#   mu <- seq(-1.5, 1.5, length.out = k)
#
#   ## Continuous block
#   X_cont <- matrix(NA_real_, nrow = n, ncol = num_signal_cont + num_noise_cont)
#   colnames(X_cont) <- c(
#     paste0("ContSig_", seq_len(num_signal_cont)),
#     paste0("ContNoise_", seq_len(num_noise_cont))
#   )
#
#   # Signal continuous vars: different means per class
#   for (j in seq_len(num_signal_cont)) {
#     X_cont[, j] <- rnorm(n, mean = mu[as.numeric(y)], sd = cont_signal_sd)
#   }
#   # Noise continuous vars: independent of y
#   for (j in seq_len(num_noise_cont)) {
#     X_cont[, num_signal_cont + j] <- rnorm(n, sd = cont_noise_sd)
#   }
#
#   ## Categorical base variable strongly associated with y
#   base_cat <- character(n)
#   for (g in seq_len(k)) {
#     idx <- which(y == levels(y)[g])
#     base_probs <- rep((1 - cat_assoc_strength) / (levels_per_cat - 1),
#                       levels_per_cat)
#     main_level <- sample(seq_len(levels_per_cat), 1)
#     base_probs[main_level] <- cat_assoc_strength
#     base_cat[idx] <- sample(
#       paste0("CatBase_L", seq_len(levels_per_cat)),
#       size = length(idx),
#       replace = TRUE,
#       prob = base_probs
#     )
#   }
#   base_cat <- factor(base_cat)
#
#   # Redundant categorical variables: noisy copies of base_cat
#   cat_list <- list(base_cat)
#   if (num_redundant_cat > 0) {
#     for (c in seq_len(num_redundant_cat)) {
#       tmp <- base_cat
#       flip_idx <- sample(seq_len(n), size = round(flip_prob_redundant * n))
#       levs <- levels(tmp)
#       for (i in flip_idx) {
#         tmp[i] <- sample(levs[levs != tmp[i]], 1)
#       }
#       cat_list[[length(cat_list) + 1]] <- droplevels(tmp)
#     }
#   }
#
#   # Pure noise categorical variables: unrelated to y or base_cat
#   if (num_cat_noise > 0) {
#     for (c in seq_len(num_cat_noise)) {
#       cat_noise <- sample(
#         paste0("CatNoise", c, "_L", seq_len(levels_per_cat)),
#         size = n,
#         replace = TRUE
#       )
#       cat_list[[length(cat_list) + 1]] <- factor(cat_noise)
#     }
#   }
#
#   X_cat <- as.data.frame(cat_list)
#   cn <- c("CatBase", paste0("CatRed_", seq_len(num_redundant_cat)),
#           paste0("CatNoise_", seq_len(num_cat_noise)))
#   colnames(X_cat) <- cn
#
#   df <- cbind(as.data.frame(X_cont), X_cat)
#   df$y <- y
#   df
# }
