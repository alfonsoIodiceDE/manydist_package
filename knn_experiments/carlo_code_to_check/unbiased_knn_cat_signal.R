## -------------------------------------------------
## Packages
## -------------------------------------------------
packages <- c(
  "mclust", "class", "manydist", "tidyverse", "tidymodels",
  "data.table", "dplyr", "MASS", "rsample", "cluster", "clusterGeneration"
)
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

## -------------------------------------------------
## Source functions
## -------------------------------------------------
source("functions_unbiased_knn.R")

## -------------------------------------------------
## Simulation settings (single place)
## -------------------------------------------------
nsample <- 10

k_true  <- 4                    # true number of classes in df$y
k_max   <- 15 

signal <- 4
noise <- 4
q <- 5

sep <- 0.05

methods <- c("g", "oh", "hl", "hla", "uind", "udep")

acc_summary <- vector("list", nsample)
ari_summary <- vector("list", nsample)

for (s in seq_len(nsample)) {
  
  ## ----- Generate data -----
  gen <- genRandomClust(numClust = k_true, sepVal = sep, numNonNoisy = signal, numNoisy = noise, numReplicate = 1)
  
  X <- gen$datList$test_1
  p <- signal + noise
  
  noise_idx  <- sort(gen$noisyList$test_1)
  signal_idx <- setdiff(seq_len(p), noise_idx)
  
  df <- as.data.frame(X)
  for (j in signal_idx) {
    df[[j]] <- as.factor(cut(X[, j], breaks = q, labels = seq_len(q)))
  }
  
  df$y <- as.factor(gen$memList$test_1)
  
  ## -------------------------------------------------
  ## 1) KNN + accuracy (supervised)
  ## -------------------------------------------------
  split_obj <- initial_split(df, prop = 0.8, strata = "y")
  tr_df     <- training(split_obj)
  ts_df     <- testing(split_obj)
  
  y_tr <- tr_df$y
  y_ts <- ts_df$y
  
  tr_x <- dplyr::select(tr_df, -y)
  ts_x <- dplyr::select(ts_df, -y)
  
  X_full <- dplyr::select(df, -y)
  
  dists_tt <- compute_all_dist_train_test(tr_x, ts_x)
  
  acc_mat <- matrix(NA_real_, nrow = k_max, ncol = length(methods),
                    dimnames = list(paste0("k=", 1:k_max), methods))
  
  for (k in 1:k_max) {
    for (m in methods) {
      pred <- knn_dist(dists_tt[[m]], y_tr, k)
      acc_mat[k, m] <- compute_accuracy(pred, y_ts)
    }
  }
  
  acc_summary[[s]] <- tibble(
    rep = s,
    method = factor(methods, levels = methods),
    accuracy = apply(acc_mat, 2, max, na.rm = TRUE)
  )
  
  ## -------------------------------------------------
  ## 2) PAM + ARI (unsupervised)
  ## -------------------------------------------------
  y_full <- df$y
  
  dists_full <- compute_all_dist_full(X_full)
  
  ari_vec <- sapply(methods, function(m) {
    pam_fit <- pam(dists_full[[m]], k = k_true, diss = TRUE)
    adjustedRandIndex(y_full, pam_fit$clustering)
  })
  
  ari_summary[[s]] <- tibble(
    rep = s,
    method = factor(names(ari_vec), levels = methods),
    ari = as.numeric(ari_vec)
  )
}

acc_summary <- bind_rows(acc_summary)
ari_summary <- bind_rows(ari_summary)

## -------------------------------------------------
## Plots
## -------------------------------------------------
# ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
#   geom_boxplot(alpha = 0.7) +
#   theme_minimal(base_size = 14) +
#   labs(x = "Distance Method", y = "Accuracy (KNN)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme(legend.position = "none")

ggplot(ari_summary, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none")

save.image("unbiased_knn_cat_signal.RData")