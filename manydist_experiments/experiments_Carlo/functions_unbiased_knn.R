compute_all_dist_full <- function(X_full) {
  list(
    g    = mdist(X_full, preset = "gower")$distance,
    oh   = mdist(X_full, preset = "euclidean_onehot")$distance,
    mod_g = mdist(X_full, preset = "mod_gower")$distance,
    dkss  = mdist(X_full, preset = "dkss")$distance,
    gudmm = mdist(X_full, preset = "gudmm")$distance,
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
    udep = mdist(X_full, preset = "unbiased_dependent")$distance,
    umix = mdist(X_full, distance_cont  = "manhattan",
                 distance_cat   = "tot_var_dist",
                 commensurable  = TRUE)$distance
  )
}

gen_mixed <- function(
    k_true,
    clustSizeEq = 50,
    numsignal = 2,
    numnoise  = 2,
    catsignal = 2,
    catnoise  = 2,
    q = 5,
    q_err = 9,
    numsep = 0.1,
    catsep = 0.5,
    seed = NULL,
    error_type = c("normal", "chisq"),
    error_df = 2,
    error_scale = 1
) {
  error_type <- match.arg(error_type)
  if (!is.null(seed)) set.seed(seed)
  
  X_num_df <- NULL
  X_cat_df <- NULL
  y_ref    <- NULL
  n        <- k_true * clustSizeEq
  
  make_balanced_y <- function(k, m) {
    factor(rep(seq_len(k), each = m), levels = seq_len(k))
  }
  
  gen_noise <- function(n, p, type = c("normal", "chisq"), df = .5, scale = 2) {
    type <- match.arg(type)
    if (p <= 0) return(matrix(numeric(0), nrow = n, ncol = 0))
    if (type == "normal") {
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
    } else {
      Z <- matrix(rchisq(n * p, df = df) - df, nrow = n, ncol = p)
    }
    Z * scale
  }
  
  discretize_to_q <- function(X, q) {
    dfc <- as.data.frame(X)
    if (ncol(dfc) == 0) return(dfc)
    
    for (j in seq_len(ncol(dfc))) {
      x <- X[, j]
      if (is.na(sd(x)) || sd(x) == 0) x <- x + rnorm(length(x), sd = 1e-8)
      cj <- cut(x, breaks = q, include.lowest = TRUE)
      
      if (all(is.na(cj))) {
        brks <- unique(quantile(x, probs = seq(0, 1, length.out = q + 1), na.rm = TRUE))
        if (length(brks) <= 2) {
          x <- x + rnorm(length(x), sd = 1e-4)
          brks <- unique(quantile(x, probs = seq(0, 1, length.out = q + 1), na.rm = TRUE))
        }
        cj <- cut(x, breaks = brks, include.lowest = TRUE)
      }
      
      cj <- droplevels(cj)
      if (nlevels(cj) < 2) {
        x <- x + rnorm(length(x), sd = 1e-4)
        brks <- unique(quantile(x, probs = seq(0, 1, length.out = q + 1), na.rm = TRUE))
        cj <- droplevels(cut(x, breaks = brks, include.lowest = TRUE))
      }
      
      dfc[[j]] <- factor(as.integer(cj), levels = seq_len(min(q, nlevels(cj))))
    }
    dfc
  }
  
  # 1) numerical signal
  if (numsignal > 0) {
    gen_num <- genRandomClust(
      numClust     = k_true,
      sepVal       = numsep,
      numNonNoisy  = numsignal,
      numNoisy     = 0,
      numReplicate = 1,
      clustszind   = 1,
      clustSizeEq  = clustSizeEq
    )
    
    X_num_sig <- gen_num$datList$test_1
    y_num     <- factor(gen_num$memList$test_1, levels = seq_len(k_true))
    
    dfn <- as.data.frame(X_num_sig)
    colnames(dfn) <- paste0("NumSig_", seq_len(ncol(dfn)))
    dfn$y <- y_num
    
    dfn <- dfn[order(dfn$y), , drop = FALSE]
    y_ref <- dfn$y
    X_num_df <- dfn
  }
  
  # 2) categorical signal
  if (catsignal > 0) {
    gen_cat <- genRandomClust(
      numClust     = k_true,
      sepVal       = catsep,
      numNonNoisy  = catsignal,
      numNoisy     = 0,
      numReplicate = 1,
      clustszind   = 1,
      clustSizeEq  = clustSizeEq
    )
    
    X_cat_sig_num <- gen_cat$datList$test_1
    y_cat         <- factor(gen_cat$memList$test_1, levels = seq_len(k_true))
    
    dfc <- discretize_to_q(X_cat_sig_num, q = q)
    colnames(dfc) <- paste0("CatSig_", seq_len(ncol(dfc)))
    dfc$y <- y_cat
    dfc <- dfc[order(dfc$y), , drop = FALSE]
    
    if (is.null(y_ref)) y_ref <- dfc$y else dfc$y <- y_ref
    X_cat_df <- dfc
  }
  
  # 3) if no signal
  if (is.null(y_ref)) {
    y_ref <- make_balanced_y(k_true, clustSizeEq)
  }
  
  # Ensure empty blocks are explicit 0-column data.frames with n rows
  if (is.null(X_num_df)) X_num_df <- data.frame(y = y_ref)[, 0, drop = FALSE]
  if (is.null(X_cat_df)) X_cat_df <- data.frame(y = y_ref)[, 0, drop = FALSE]
  
  # 4) numerical noise
  if (numnoise > 0) {
    Xn <- gen_noise(n = n, p = numnoise, type = error_type, df = error_df, scale = error_scale)
    dn <- as.data.frame(Xn)
    colnames(dn) <- paste0("NumNoise_", seq_len(ncol(dn)))
    dn <- dn[order(y_ref), , drop = FALSE]
    X_num_df <- cbind(X_num_df, dn)
  }
  
  # -------------------------------------------------
  # 5) Add CATEGORICAL NOISE via (numeric noise -> discretize to q)
  #     - same error generation as numerical part
  #     - then discretization to q levels
  # -------------------------------------------------
  if (catnoise > 0) {
    # 1) generate numeric noise with same mechanism as numerical noise
    Xc_num <- gen_noise(
      n = n, p = catnoise,
      type = error_type, df = error_df, scale = error_scale
    )
    
    # 2) discretize to q-level factors
    dc <- discretize_to_q(Xc_num, q = q_err)
    
    # 3) name columns and align with y_ref ordering
    colnames(dc) <- paste0("CatNoise_", seq_len(ncol(dc)))
    
    # attach y, order, then drop y (to mimic your other blocks)
    dc$y <- y_ref
    dc <- dc[order(dc$y), , drop = FALSE]
    dc <- dc[, setdiff(names(dc), "y"), drop = FALSE]
    
    # 4) append to existing categorical df (or initialize if needed)
    if (is.null(X_cat_df)) {
      X_cat_df <- dc
    } else {
      X_cat_df <- cbind(X_cat_df, dc)
    }
  }
  
  # 6) build df
  df <- data.frame(y = y_ref)
  
  if (ncol(X_num_df) > 0) df <- cbind(df, X_num_df)
  if (ncol(X_cat_df) > 0) df <- cbind(df, X_cat_df)
  
  df <- df[, c(setdiff(names(df), "y"), "y"), drop = FALSE]
  
  num_cols <- grep("^Num", names(df), value = TRUE)
  cat_cols <- grep("^Cat", names(df), value = TRUE)
  
  return(list(
    df       = df,
    X_num    = X_num_df,
    X_cat    = X_cat_df,
    y        = y_ref,
    num_cols = num_cols,
    cat_cols = cat_cols
  ))
}