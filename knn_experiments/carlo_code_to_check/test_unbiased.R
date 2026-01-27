## -------------------------------------------------
## Packages
## -------------------------------------------------
packages <- c("mclust", "class", "tidyverse", "tidymodels", "data.table", "dplyr", "MASS", "rsample", "cluster", "clusterGeneration", "parallelDist","klaR","arules","FD","StatMatch","clustMixType")

to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

# options(repos = c(CRAN = "https://cloud.r-project.org"))
# remotes::install_github("alfonsoIodiceDE/manydist_package", subdir = "manydist")
library(manydist)

## -------------------------------------------------
## Source functions
## -------------------------------------------------
source("../knn_experiments/carlo_code_to_check/functions_unbiased_knn.R")

## -------------------------------------------------
## Simulation settings (single place)
## -------------------------------------------------
nsample <- 10

k_true  <- 4   # true number of classes in df$y
k_max   <- 15  # max value for k for knn

#methods_knn <- c("g", "oh", "hl", "hla", "uind", "udep") # fix manydist before running knn
methods_pam <- c("g", "oh", "hl", "hla", "uind", "udep") #add the other methods

acc_summary <- vector("list", nsample)
ari_summary_num <- vector("list", nsample)
ari_summary_cat <- vector("list", nsample)
ari_summary     <- vector("list", nsample)

# Mixed data generation controls
numsignal <- 3
numnoise  <- 3
catsignal <- 3
catnoise  <- 3
q <- 9

numsep <- 0.1
catsep <- 0.6

# Fixed cluster size (so n = k_true * clustSizeEq)
clustSizeEq <- 50

## -------------------------------------------------
## Main loop
## -------------------------------------------------
for (s in seq_len(nsample)) {

  ## ----- Generate mixed data (single wrapper call) -----
  out <- gen_mixed(
    k_true      = k_true,
    clustSizeEq = clustSizeEq,
    numsignal   = numsignal,
    numnoise    = numnoise,
    catsignal   = catsignal,
    catnoise    = catnoise,
    q           = q,
    numsep      = numsep,
    catsep      = catsep,
    seed        = s
  )

  df  <- out$df
  dfn <- out$X_num   # may be NULL if no numeric cols
  dfc <- out$X_cat   # may be NULL if no categorical cols

  # ## -------------------------------------------------
  # ## 1) KNN + accuracy (supervised)
  # ## -------------------------------------------------
  #
  # ## --- Numerical-only accuracy (if numeric block exists) ---
  # if (!is.null(dfn) && ncol(dfn) > 1) {
  #   split_obj <- initial_split(dfn, prop = 0.8, strata = "y")
  #   tr_df     <- training(split_obj)
  #   ts_df     <- testing(split_obj)
  #
  #   y_tr <- tr_df$y
  #   y_ts <- ts_df$y
  #
  #   tr_x <- tr_df[, setdiff(names(tr_df), "y"), drop = FALSE]
  #   ts_x <- ts_df[, setdiff(names(ts_df), "y"), drop = FALSE]
  #
  #   dists_tt_num <- compute_all_dist_train_test(tr_x, ts_x)
  #
  #   acc_mat_num <- matrix(
  #     NA_real_,
  #     nrow = k_max,
  #     ncol = length(methods_knn),
  #     dimnames = list(paste0("k=", 1:k_max), methods_knn)
  #   )
  #
  #   for (k in 1:k_max) {
  #     for (m in methods_knn) {
  #       pred <- knn_dist(dists_tt_num[[m]], y_tr, k)
  #       acc_mat_num[k, m] <- compute_accuracy(pred, y_ts)
  #     }
  #   }
  #
  #   acc_summary_num[[s]] <- tibble(
  #     rep       = s,
  #     method_knn = factor(methods_knn, levels = methods_knn),
  #     accuracy  = apply(acc_mat_num, 2, max, na.rm = TRUE)
  #   )
  #
  # } else {
  #   acc_summary_num[[s]] <- tibble(
  #     rep       = s,
  #     method_knn = factor(methods_knn, levels = methods_knn),
  #     accuracy  = NA_real_
  #   )
  # }
  #
  # ## --- Categorical-only accuracy (if categorical block exists) ---
  # if (!is.null(dfc) && ncol(dfc) > 1) {
  #   split_obj <- initial_split(dfc, prop = 0.8, strata = "y")
  #   tr_df     <- training(split_obj)
  #   ts_df     <- testing(split_obj)
  #
  #   y_tr <- tr_df$y
  #   y_ts <- ts_df$y
  #
  #   tr_x <- tr_df[, setdiff(names(tr_df), "y"), drop = FALSE]
  #   ts_x <- ts_df[, setdiff(names(ts_df), "y"), drop = FALSE]
  #
  #   dists_tt_cat <- compute_all_dist_train_test(tr_x, ts_x)
  #
  #   acc_mat_cat <- matrix(
  #     NA_real_,
  #     nrow = k_max,
  #     ncol = length(methods_knn),
  #     dimnames = list(paste0("k=", 1:k_max), methods_knn)
  #   )
  #
  #   for (k in 1:k_max) {
  #     for (m in methods_knn) {
  #       pred <- knn_dist(dists_tt_cat[[m]], y_tr, k)
  #       acc_mat_cat[k, m] <- compute_accuracy(pred, y_ts)
  #     }
  #   }
  #
  #   acc_summary_cat[[s]] <- tibble(
  #     rep       = s,
  #     method_knn = factor(methods_knn, levels = methods_knn),
  #     accuracy  = apply(acc_mat_cat, 2, max, na.rm = TRUE)
  #   )
  #
  # } else {
  #   acc_summary_cat[[s]] <- tibble(
  #     rep       = s,
  #     method_knn = factor(methods_knn, levels = methods_knn),
  #     accuracy  = NA_real_
  #   )
  # }
  #
  # ## --- Mixed accuracy (always, from df) ---
  # split_obj <- initial_split(df, prop = 0.8, strata = "y")
  # tr_df     <- training(split_obj)
  # ts_df     <- testing(split_obj)
  #
  # y_tr <- tr_df$y
  # y_ts <- ts_df$y
  #
  # tr_x <- tr_df[, setdiff(names(tr_df), "y"), drop = FALSE]
  # ts_x <- ts_df[, setdiff(names(ts_df), "y"), drop = FALSE]
  #
  # dists_tt_mix <- compute_all_dist_train_test(tr_x, ts_x)
  #
  # acc_mat <- matrix(
  #   NA_real_,
  #   nrow = k_max,
  #   ncol = length(methods_knn),
  #   dimnames = list(paste0("k=", 1:k_max), methods_knn)
  # )
  #
  # for (k in 1:k_max) {
  #   for (m in methods_knn) {
  #     pred <- knn_dist(dists_tt_mix[[m]], y_tr, k)
  #     acc_mat[k, m] <- compute_accuracy(pred, y_ts)
  #   }
  # }
  #
  # acc_summary[[s]] <- tibble(
  #   rep       = s,
  #   method_knn = factor(methods_knn, levels = methods_knn),
  #   accuracy  = apply(acc_mat, 2, max, na.rm = TRUE)
  # )

  ## -------------------------------------------------
  ## 2) PAM + ARI (unsupervised)
  ## -------------------------------------------------

  ## --- Numerical-only ARI (if numeric block exists) ---
  if (!is.null(dfn) && ncol(dfn) > 1) {
    X_num <- dfn[, setdiff(names(dfn), "y"), drop = FALSE]
    y_num <- dfn$y

    dists_full_num <- compute_all_dist_full(X_num)

    ari_vec_num <- sapply(methods_pam, function(m) {
      pam_fit <- pam(dists_full_num[[m]], k = k_true, diss = TRUE)
      adjustedRandIndex(y_num, pam_fit$clustering)
    })

    ari_summary_num[[s]] <- tibble(
      rep    = s,
      method_pam = factor(names(ari_vec_num), levels = methods_pam),
      ari    = as.numeric(ari_vec_num)
    )
  } else {
    ari_summary_num[[s]] <- tibble(
      rep    = s,
      method_pam = factor(methods_pam, levels = methods_pam),
      ari    = NA_real_
    )
  }

  ## --- Categorical-only ARI (if categorical block exists) ---
  if (!is.null(dfc) && ncol(dfc) > 1) {
    X_cat <- dfc[, setdiff(names(dfc), "y"), drop = FALSE]
    y_cat <- dfc$y

    dists_full_cat <- compute_all_dist_full(X_cat)

    ari_vec_cat <- sapply(methods_pam, function(m) {
      pam_fit <- pam(dists_full_cat[[m]], k = k_true, diss = TRUE)
      adjustedRandIndex(y_cat, pam_fit$clustering)
    })

    ari_summary_cat[[s]] <- tibble(
      rep    = s,
      method_pam = factor(names(ari_vec_cat), levels = methods_pam),
      ari    = as.numeric(ari_vec_cat)
    )
  } else {
    ari_summary_cat[[s]] <- tibble(
      rep    = s,
      method_pam = factor(methods_pam, levels = methods_pam),
      ari    = NA_real_
    )
  }

  ## --- Mixed ARI (always, from df) ---
  y_full <- df$y
  X_full <- df[, setdiff(names(df), "y"), drop = FALSE]

  dists_full_mix <- compute_all_dist_full(X_full)

  ari_vec <- sapply(methods_pam, function(m) {
    pam_fit <- pam(dists_full_mix[[m]], k = k_true, diss = TRUE)
    adjustedRandIndex(y_full, pam_fit$clustering)
  })

  ari_summary[[s]] <- tibble(
    rep    = s,
    method_pam = factor(names(ari_vec), levels = methods_pam),
    ari    = as.numeric(ari_vec)
  )
}

## -------------------------------------------------
## Bind results
## -------------------------------------------------
acc_summary      <- bind_rows(acc_summary)
ari_summary_num  <- bind_rows(ari_summary_num)
ari_summary_cat  <- bind_rows(ari_summary_cat)
ari_summary      <- bind_rows(ari_summary)

## -------------------------------------------------
## Plots
## -------------------------------------------------
# accplot <- acc_summary_num
# variant <- paste(numsignal, "Numerical signal +", numnoise, "num. noise", sep = " ")
#
# ggplot(accplot, aes(x = method_knn, y = accuracy, fill = method_knn)) +
#   geom_boxplot(alpha = 0.7) +
#   theme_minimal(base_size = 14) +
#   labs(x = "Distance Method", y = "Accuracy (KNN)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme(legend.position = "none") +
#   ggtitle(variant)
#
# accplot <- acc_summary_cat
# variant <- paste(catsignal, "Categorical signal +", catnoise, "cat. noise", sep = " ")
#
# ggplot(accplot, aes(x = method_knn, y = accuracy, fill = method_knn)) +
#   geom_boxplot(alpha = 0.7) +
#   theme_minimal(base_size = 14) +
#   labs(x = "Distance Method", y = "Accuracy (KNN)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme(legend.position = "none") +
#   ggtitle(variant)
#
# accplot <- acc_summary
# variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
#                  catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")
#
# ggplot(accplot, aes(x = method_knn, y = accuracy, fill = method_knn)) +
#   geom_boxplot(alpha = 0.7) +
#   theme_minimal(base_size = 14) +
#   labs(x = "Distance Method", y = "Accuracy (KNN)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme(legend.position = "none") +
#   ggtitle(variant)

ariplot <- ari_summary_num
variant <- paste(numsignal, "Numerical signal +", numnoise, "num. noise", sep = " ")

ggplot(ariplot, aes(x = method_pam, y = ari, fill = method_pam)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

ariplot <- ari_summary_cat
variant <- paste(catsignal, "Categorical signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method_pam, y = ari, fill = method_pam)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method_pam, y = ari, fill = method_pam)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)
