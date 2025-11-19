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
