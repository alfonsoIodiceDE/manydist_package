# needs: purrr
delta_knn_ba <- function(
    D, labels,
    pi_nn   = 0.1,                      # neighbor proportion in (0,1]
    decision = c("posterior", "prior_corrected"),
    nn_order = NULL                     # optional precomputed n x (n-1) neighbor order
) {
  decision <- match.arg(decision)
  stopifnot(is.matrix(D), nrow(D) == ncol(D))
  n <- nrow(D); stopifnot(length(labels) == n)
  labels <- droplevels(factor(labels))
  L <- levels(labels); G <- length(L)
  
  # Early return if <2 classes
  Delta <- matrix(0, G, G, dimnames = list(L, L))
  if (G < 2) return(Delta)
  
  # --- Use external nn_order if provided, else compute ---
  if (!is.null(nn_order)) {
    # basic checks
    if (!is.matrix(nn_order) || nrow(nn_order) != n || ncol(nn_order) != (n - 1))
      stop("nn_order must be an n x (n-1) integer matrix matching D.")
    # enforce integers within 1..n and no self indices per row
    if (any(nn_order < 1 | nn_order > n, na.rm = TRUE))
      stop("nn_order contains indices outside [1, n].")
    # remove any accidental self indices (just in case)
    for (i in seq_len(n)) {
      nn_order[i, ] <- nn_order[i, nn_order[i, ] != i]
      # pad if any removed (rare); recompute from D for that row
      if (length(nn_order[i, ]) < (n - 1)) {
        ord <- order(D[, i], na.last = NA)
        nn_order[i, ] <- ord[ord != i][seq_len(n - 1)]
      }
    }
    NN_order <- nn_order
  } else {
    NN_order <- matrix(NA_integer_, n, n - 1)
    for (i in seq_len(n)) {
      ord <- order(D[, i], na.last = NA)
      ord <- ord[ord != i]
      NN_order[i, ] <- ord
    }
  }
  
  # helpers ------------------------------------------------------
  k_from_prop <- function(m) {
    # m = |{a,b}|, neighbors exclude self -> m-1
    k <- floor(pi_nn * (m - 1L))
    max(1L, min(m - 1L, as.integer(k)))
  }
  
  # mean of "is target" over first k neighbors restricted to 'allowed' set
  knn_prop_target <- function(i, allowed_mask, k, target_level) {
    row_ord <- NN_order[i, ]
    # filter to allowed candidates
    cand <- row_ord[allowed_mask[row_ord]]
    if (!length(cand)) return(NA_real_)
    k_use <- min(k, length(cand))
    idx <- cand[seq_len(k_use)]
    mean(labels[idx] == target_level)
  }
  
  # --------------------------------------------------------------
  for (a in seq_len(G - 1L)) {
    for (b in seq.int(a + 1L, G)) {
      pair_lvls <- c(L[a], L[b])
      pair_idx  <- which(labels %in% pair_lvls)
      m <- length(pair_idx)
      if (m <= 1L) next
      
      k_use <- k_from_prop(m)
      
      # neighbors restricted to {a,b}
      allowed <- rep(FALSE, n); allowed[pair_idx] <- TRUE
      
      # soft scores: p_a(i) for each i in the pair
      p_a <- purrr::map_dbl(pair_idx, ~ knn_prop_target(.x, allowed, k_use, L[a]))
      
      # decision threshold
      thr <- if (decision == "posterior") 0.5 else mean(labels[pair_idx] == L[a])
      
      # hard predictions (ties go to 'a' via >=)
      yhat_is_a <- (p_a >= thr)
      
      # true labels
      y_is_a <- (labels[pair_idx] == L[a])
      
      # recalls and BA
      rec_a <- if (any(y_is_a)) mean(yhat_is_a[y_is_a], na.rm = TRUE) else 0
      rec_b <- if (any(!y_is_a)) mean(!yhat_is_a[!y_is_a], na.rm = TRUE) else 0
      BA    <- 0.5 * (rec_a + rec_b)
      
      # map BA ∈ [0.5,1] to Δ ∈ [0,1]
      val <- max(0, min(1, 2 * (BA - 0.5)))
      
      Delta[a, b] <- val
      Delta[b, a] <- val
    }
  }
  
  diag(Delta) <- 0
  Delta
}
