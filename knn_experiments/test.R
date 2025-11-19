packages <- c("mclust", "class", "manydist", "tidyverse", "tidymodels", "data.table", "dplyr","MASS", "rsample")
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

n <- 100
num_continuous <- 5
num_noise <- 5
noise_strength <- 10
num_categorical <- 5
num_categorical_noise <- 5
target_corr <- 0.7
cat_assoc_strength <- 0.4

y <- factor(sample(c("A", "B", "C"), size = n, replace = TRUE))

df <- genmixdata(n = n,
                 num_continuous = num_continuous, 
                 target_corr =  target_corr, 
                 y = y, 
                 num_noise = num_noise, 
                 noise_strength = noise_strength,
                 num_categorical = num_categorical, 
                 levels_per_cat = 3, 
                 cat_assoc_strength = cat_assoc_strength,
                 num_categorical_noise = num_categorical_noise)

acc_summary <- data.frame()
nsample <- 10

for (s in  1:nsample) {
  df_split <- initial_split(df, prop = 0.7)
  tr_df <- training(df_split)
  ts_df <- testing(df_split)
  
  y_tr <- y[df_split$in_id]
  y_ts <- y[-df_split$in_id]
  
  tr_num <- tr_df[, 1:(num_continuous + num_noise)]
  ts_num <- ts_df[, 1:(num_continuous + num_noise)]
  tr_cat <- tr_df[, (num_continuous + num_noise + 1):(num_continuous + num_noise + num_categorical)]
  ts_cat <- ts_df[, (num_continuous + num_noise + 1):(num_continuous + num_noise + num_categorical)]
  
  # Numerical distances
  dist_mana <- mana(tr_num, ts_num, sig = diag(ncol(tr_num)))
  sigma_sv <- diag(sapply(1:ncol(tr_num), function(j) dann_metric_sv(tr_num[, j], y_tr)))
  dist_lda_sv_add <- mana(tr_num, ts_num, sig = sigma_sv)
  sigma <- dann_metric(tr_num, y_tr)
  dist_lda <- maha(tr_num, ts_num, sig = sigma)
  dist_lda_sv <- maha(tr_num, ts_num, sig = sigma_sv)
  
  # Categorical distances
  stvdf <- cdist(x = tr_cat, validate_x = ts_cat, response = y_tr, method = "supervised_full", commensurable = TRUE)$distance_mat
  stvd <- cdist(x = tr_cat, validate_x = ts_cat, response = y_tr, method = "supervised", commensurable = TRUE)$distance_mat
  #match <- cdist(x = tr_cat, validate_x = ts_cat, method = "euc", commensurable = FALSE)$distance_mat
  
  # Combine distances
  # dist_sup <- (ncol(tr_df)/ncol(tr_num))*dist_lda_sv + (ncol(tr_df)/ncol(tr_cat))*t(stvd)
  # dist_sup_add <- (ncol(tr_df)/ncol(tr_num))*dist_lda_sv_add + (ncol(tr_df)/ncol(tr_cat))*t(stvd)
  # dist_supf <- (ncol(tr_df)/ncol(tr_num))*dist_lda + (ncol(tr_df)/ncol(tr_cat))*t(stvdf)
  # dist_nsup <- t(cran_mdist(tr_df, ts_df, preset = "euclidean_onehot"))
  # dist_g <- t(cran_mdist(tr_df, ts_df, preset = "gower"))
  
  dist_sup <- dist_lda_sv + t(stvd)
  dist_sup_add <- dist_lda_sv_add + t(stvd)
  dist_supf <- dist_lda + t(stvdf)
  dist_nsup <- t(cran_mdist(tr_df, ts_df, preset = "euclidean_onehot"))
  dist_g <- t(cran_mdist(tr_df, ts_df, preset = "gower"))
  
  # Accuracy for each k
  k_max <- 15
  acc <- matrix(0, k_max, 5)
  for (k in 1:k_max) {
    acc[k, 1] <- compute_accuracy(y_ts, knn_dist(dist_sup, y_tr, k))
    acc[k, 2] <- compute_accuracy(y_ts, knn_dist(dist_sup_add, y_tr, k))
    acc[k, 3] <- compute_accuracy(y_ts, knn_dist(dist_supf, y_tr, k))
    acc[k, 4] <- compute_accuracy(y_ts, knn_dist(dist_nsup, y_tr, k))
    acc[k, 5] <- compute_accuracy(y_ts, knn_dist(dist_g, y_tr, k))
  }
  
  # Store best accuracy for each method
  acc_summary <- rbind(acc_summary, data.frame(
    rep = s, 
    method = c("sup", "sup_add", "supf", "naive", "gower"),
    accuracy = apply(acc, 2, max)
  ))
}

ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Best k-NN Accuracy",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme(legend.position = "none")
