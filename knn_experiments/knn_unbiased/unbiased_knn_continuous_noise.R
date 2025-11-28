packages <- c("mclust", "class", "manydist", "tidyverse", "tidymodels", "data.table", "dplyr","MASS", "rsample")
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

source("functions_unbiased_knn.R")
n <- 1000
num_continuous <- 10
num_noise <- 10
num_categorical <- 5
num_categorical_noise <- 10

# continuous:
target_corr <- rep(0.7, num_continuous) 

# categorical: 
cat_assoc_strength <- rep(0.2, num_categorical)

y <- factor(sample(c("A", "B", "C"), size = n, replace = TRUE))

df <- genmixdata(
  n = n,
  num_continuous = num_continuous,
  target_corr = target_corr,
  y = y,
  num_noise = num_noise,
  noise_strength = 5,
  num_categorical = num_categorical,
  levels_per_cat = 3,
  cat_assoc_strength = cat_assoc_strength,
  cat_redundancy = TRUE,            
  cat_redundancy_noise_prob = 0.1,   
  num_categorical_noise = num_categorical_noise
)

df$y <- y 

acc_summary <- data.frame()
nsample <- 100
k_max <- 15

for (s in 1:nsample) {
  df$y <- y 
  df_split <- initial_split(df, prop = 0.7, strata = y)
  df <- df[, !(names(df) == "y")]
  tr_df <- training(df_split)
  ts_df <- testing(df_split)
  
  y_tr <- tr_df$y
  y_ts <- ts_df$y
  
  tr_df$y <- NULL
  ts_df$y <- NULL
  
  dist_g   <- t(mdist(tr_df, validate_x = ts_df, preset = "gower"))
  dist_oh  <- t(mdist(tr_df, validate_x = ts_df, preset = "euclidean_onehot"))
  dist_udep <- t(mdist(tr_df, validate_x = ts_df, preset = "unbiased_dependent"))
  dist_uind <- t(mdist(tr_df, validate_x = ts_df, distance_cont = "manhattan", distance_cat  = "cat_dis", commensurable = TRUE))
  
  acc <- matrix(0, k_max, 4)
  for (k in 1:k_max) {
    acc[k, 1] <- compute_accuracy(y_ts, knn_dist(dist_g,   y_tr, k))
    acc[k, 2] <- compute_accuracy(y_ts, knn_dist(dist_oh,  y_tr, k))
    acc[k, 3] <- compute_accuracy(y_ts, knn_dist(dist_udep, y_tr, k))
    acc[k, 4] <- compute_accuracy(y_ts, knn_dist(dist_uind, y_tr, k))
  }
  
  acc_summary <- rbind(acc_summary, data.frame(
    rep = s,
    method = c("g", "oh", "udep", "uind"),
    accuracy = apply(acc, 2, max) 
  ))
}

ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none")

save.image("unbiased_knn_continuous_noise.RData")