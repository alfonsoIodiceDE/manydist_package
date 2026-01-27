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
compute_all_dist_train_test <- function(tr_df, ts_df) {
  list(
    g    = t(mdist(tr_df, new_data = ts_df, preset = "gower")$distance |> as.matrix()),
    #oh   = t(mdist(tr_df, validate_x = ts_df, preset = "euclidean_onehot")),
    hl   = t(mdist(tr_df, new_data = ts_df, preset = "custom",
                   distance_cont = "euclidean", distance_cat = "HLeucl",
                   scaling_cont  = "std")$distance |> as.matrix()),
    hla  = t(mdist(tr_df, new_data = ts_df, preset = "custom",
                   distance_cont = "manhattan", distance_cat = "HL",
                   scaling_cont  = "std")$distance |> as.matrix()),
    uind = t(mdist(tr_df, new_data = ts_df,
                   distance_cont = "manhattan", distance_cat = "cat_dis",
                   commensurable = TRUE)$distance |> as.matrix()),
    udep = t(mdist(tr_df, new_data = ts_df, preset = "unbiased_dependent")$distance |> as.matrix())
  )
}

compute_all_dist_full <- function(X_full) {
  list(
    g    = mdist(X_full, preset = "gower")$distance,
    oh   = mdist(X_full, preset = "euclidean_onehot")$distance,
    hl   = mdist(X_full, preset = "custom",
                           distance_cont = "euclidean",
                           distance_cat  = "HLeucl",
                           scaling_cont  = "std")$distance,
    hla  = mdist(X_full, preset = "custom",
                           distance_cont = "manhattan",
                           distance_cat  = "HL",
                           scaling_cont  = "std")$distance,
    uind = mdist(X_full,
                           distance_cont  = "manhattan",
                           distance_cat   = "cat_dis",
                           commensurable  = TRUE)$distance,
    udep = mdist(X_full, preset = "unbiased_dependent")$distance
  )
}

gen_mixed <- function(
    k_true,
    clustSizeEq = 50,
    numsignal = 2,
    numnoise  = 2,
    catsignal = 2,
    catnoise  = 2,
    q = 9,
    numsep = 0.1,
    catsep = 0.5,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # helper to make balanced y when needed
  make_balanced_y <- function(k, m) {
    factor(rep(seq_len(k), each = m), levels = seq_len(k))
  }

  # helper: discretize numeric matrix to q-level factors
  discretize_to_q <- function(X, q) {
    dfc <- as.data.frame(X)
    for (j in seq_len(ncol(dfc))) {
      dfc[[j]] <- factor(cut(X[, j], breaks = q, labels = seq_len(q)))
    }
    dfc
  }

  # -------------------------------------------------
  # 0) Determine n and y_ref
  # -------------------------------------------------
  y_ref <- NULL
  n <- k_true * clustSizeEq

  # -------------------------------------------------
  # 1) Numerical block (use genRandomClust if possible)
  # -------------------------------------------------
  X_num_df <- NULL

  if (numsignal > 0) {
    gen_num <- genRandomClust(
      numClust     = k_true,
      sepVal       = numsep,
      numNonNoisy  = numsignal,
      numNoisy     = numnoise,
      numReplicate = 1,
      clustszind   = 1,
      clustSizeEq  = clustSizeEq
    )

    X_num <- gen_num$datList$test_1
    y_num <- factor(gen_num$memList$test_1, levels = seq_len(k_true))

    dfn <- as.data.frame(X_num)
    colnames(dfn) <- paste0("Num_", seq_len(ncol(dfn)))
    dfn$y <- y_num

    # sort by cluster
    dfn <- dfn[order(dfn$y), , drop = FALSE]

    y_ref <- dfn$y
    X_num_df <- dfn
  } else {
    # no numerical signal: create y_ref later if needed;
    # if numnoise>0 generate pure noise once y_ref is known
    X_num_df <- NULL
  }

  # -------------------------------------------------
  # 2) Categorical block (use genRandomClust if possible)
  # -------------------------------------------------
  X_cat_df <- NULL

  if (catsignal > 0) {
    gen_cat <- genRandomClust(
      numClust     = k_true,
      sepVal       = catsep,
      numNonNoisy  = catsignal,
      numNoisy     = catnoise,
      numReplicate = 1,
      clustszind   = 1,
      clustSizeEq  = clustSizeEq
    )

    X_cat <- gen_cat$datList$test_1
    y_cat <- factor(gen_cat$memList$test_1, levels = seq_len(k_true))

    dfc <- discretize_to_q(X_cat, q = q)
    colnames(dfc) <- paste0("Cat_", seq_len(ncol(dfc)))
    dfc$y <- y_cat

    # sort by cluster
    dfc <- dfc[order(dfc$y), , drop = FALSE]

    if (is.null(y_ref)) {
      # if numerical block didn't define y, use the categorical y
      y_ref <- dfc$y
    } else {
      # enforce same y across blocks
      dfc$y <- y_ref
    }

    X_cat_df <- dfc
  } else {
    # no categorical signal: create pure noise factors later if catnoise>0
    X_cat_df <- NULL
  }

  # -------------------------------------------------
  # 3) If still no y_ref, create balanced y
  #    (this happens only if numsignal==0 AND catsignal==0)
  # -------------------------------------------------
  if (is.null(y_ref)) {
    y_ref <- make_balanced_y(k_true, clustSizeEq)
  }

  # -------------------------------------------------
  # 4) If numsignal==0, generate independent numerical noise (if requested)
  # -------------------------------------------------
  if (numsignal == 0 && numnoise > 0) {
    Xn <- matrix(rnorm(n * numnoise), nrow = n, ncol = numnoise)
    dfn <- as.data.frame(Xn)
    colnames(dfn) <- paste0("NumNoise_", seq_len(ncol(dfn)))
    dfn$y <- y_ref
    dfn <- dfn[order(dfn$y), , drop = FALSE]
    X_num_df <- dfn
  }

  # -------------------------------------------------
  # 5) If catsignal==0, generate independent categorical noise (if requested)
  # -------------------------------------------------
  if (catsignal == 0 && catnoise > 0) {
    levs <- as.character(seq_len(q))
    dfc <- as.data.frame(
      replicate(catnoise, factor(sample(levs, n, replace = TRUE), levels = levs),
                simplify = FALSE)
    )
    colnames(dfc) <- paste0("CatNoise_", seq_len(ncol(dfc)))
    dfc$y <- y_ref
    dfc <- dfc[order(dfc$y), , drop = FALSE]
    X_cat_df <- dfc
  }

  # -------------------------------------------------
  # 6) Build the final mixed data.frame
  # -------------------------------------------------
  # drop y from blocks and cbind
  X_num_only <- if (!is.null(X_num_df)) X_num_df[, setdiff(names(X_num_df), "y"), drop = FALSE] else NULL
  X_cat_only <- if (!is.null(X_cat_df)) X_cat_df[, setdiff(names(X_cat_df), "y"), drop = FALSE] else NULL

  X_full <- cbind.data.frame(X_num_only, X_cat_only)
  df <- cbind.data.frame(X_full, y = y_ref)

  return(list(
    df    = df,
    X_num = X_num_df,
    X_cat = X_cat_df,
    y     = y_ref
  ))
}
