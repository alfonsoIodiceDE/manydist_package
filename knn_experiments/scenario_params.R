## scenario_params.R
## Scenarios aligned with:
## - unbiased_knn_continuous.R
## - unbiased_knn_continuous_noise.R
## - unbiased_knn_categorical.R
## - unbiased_knn_categorical_noise.R
## - unbiased_knn_balanced.R
## - unbiased_knn_balanced_noise.R

scenarios <- list(
  # unbiased_knn_continuous.R
  continuous = list(
    n                    = 1000,
    num_continuous       = 10,
    target_corr          = 0.7,
    num_noise            = 5,
    noise_strength       = 1,
    num_categorical      = 5,
    levels_per_cat       = 5,
    cat_assoc_strength   = 0.2,
    cat_redundancy       = TRUE,
    cat_redundancy_noise_prob = 0.1,
    num_categorical_noise = 5,
    n_classes            = 3
  ),

  # unbiased_knn_continuous_noise.R
  continuous_noise = list(
    n                    = 1000,
    num_continuous       = 10,
    target_corr          = 0.7,
    num_noise            = 10,
    noise_strength       = 5,
    num_categorical      = 5,
    levels_per_cat       = 5,
    cat_assoc_strength   = 0.2,
    cat_redundancy       = TRUE,
    cat_redundancy_noise_prob = 0.1,
    num_categorical_noise = 10,
    n_classes            = 3
  ),

  # unbiased_knn_categorical.R
  categorical = list(
    n                    = 1000,
    num_continuous       = 5,
    target_corr          = 0.2,
    num_noise            = 5,
    noise_strength       = 1,
    num_categorical      = 10,
    levels_per_cat       = 5,
    cat_assoc_strength   = 0.7,
    cat_redundancy       = TRUE,
    cat_redundancy_noise_prob = 0.1,
    num_categorical_noise = 5,
    n_classes            = 3
  ),

  # unbiased_knn_categorical_noise.R
  categorical_noise = list(
    n                    = 1000,
    num_continuous       = 5,
    target_corr          = 0.2,
    num_noise            = 10,
    noise_strength       = 5,
    num_categorical      = 10,
    levels_per_cat       = 5,
    cat_assoc_strength   = 0.7,
    cat_redundancy       = TRUE,
    cat_redundancy_noise_prob = 0.1,
    num_categorical_noise = 10,
    n_classes            = 3
  ),

  # unbiased_knn_balanced.R
  balanced = list(
    n                    = 1000,
    num_continuous       = 5,
    target_corr          = 0.7,
    num_noise            = 5,
    noise_strength       = 1,
    num_categorical      = 5,
    levels_per_cat       = 5,
    cat_assoc_strength   = 0.7,
    cat_redundancy       = TRUE,
    cat_redundancy_noise_prob = 0.1,
    num_categorical_noise = 5,
    n_classes            = 3
  ),

  # unbiased_knn_balanced_noise.R
  balanced_noise = list(
    n                    = 1000,
    num_continuous       = 5,
    target_corr          = 0.7,
    num_noise            = 10,
    noise_strength       = 5,
    num_categorical      = 5,
    levels_per_cat       = 5,
    cat_assoc_strength   = 0.7,
    cat_redundancy       = TRUE,
    cat_redundancy_noise_prob = 0.1,
    num_categorical_noise = 10,
    n_classes            = 3
  )
)
