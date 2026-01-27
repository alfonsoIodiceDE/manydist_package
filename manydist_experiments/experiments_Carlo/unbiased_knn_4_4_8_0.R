nsample <- 50
k_true  <- 4   # true number of classes in df$y
numsignal <- 4
numnoise  <- 4
catsignal <- 8
catnoise  <- 0
q <- 9
numsep <- 0.1
catsep <- 0.6
clustSizeEq <- 50 # fixed cluster size (so n = k_true * clustSizeEq)

methods <- c("g", "oh", "mod_g", "dkss", "gudmm", "hl", "hla", "uind", "udep", "umix") 
ari_summary_num <- vector("list", nsample)
ari_summary_cat <- vector("list", nsample)
ari_summary     <- vector("list", nsample)

## main loop
for (s in seq_len(nsample)) {
  
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
    seed        = s,
    error_type = "normal"
  )
  
  df  <- out$df
  dfn <- out$X_num   # may be NULL if no numeric cols
  dfc <- out$X_cat   # may be NULL if no categorical cols
  
  y_full <- df$y
  X_full <- df[, setdiff(names(df), "y"), drop = FALSE]
  
  dists_full_mix <- compute_all_dist_full(X_full)
  
  ari_vec <- sapply(methods, function(m) {
    pam_fit <- pam(dists_full_mix[[m]], k = k_true, diss = TRUE)
    adjustedRandIndex(y_full, pam_fit$clustering)
  })
  
  ari_summary[[s]] <- tibble(
    rep    = s,
    method = factor(names(ari_vec), levels = methods),
    ari    = as.numeric(ari_vec)
  )
}

ari_summary      <- bind_rows(ari_summary)

## plot
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)


save.image("unbiased_knn_4_4_8_0.RData")