# needs: purrr
delta_int_knn <- function(
    D, labels,
    pi_nn   = 0.1,                      # neighbor proportion in (0,1]
    decision = c("posterior", "prior_corrected"),
    nn_order = NULL,                    # optional precomputed n x (n-1) neighbor order
    score = c("ba", "logloss"),#, "brier" NEW: how to score alignment
    eps = 1e-6                          # NEW: prob clipping for logloss
) {
  decision <- match.arg(decision)
  score <- match.arg(score)
  stopifnot(is.matrix(D), nrow(D) == ncol(D))
  n <- nrow(D); stopifnot(length(labels) == n)
  labels <- droplevels(factor(labels))
  L <- levels(labels); G <- length(L)

  # Early return if <2 classes
  Delta <- matrix(0, G, G, dimnames = list(L, L))
  if (G < 2) return(Delta)

  # --- Use external nn_order if provided, else compute ---
  if (!is.null(nn_order)) {
    if (!is.matrix(nn_order) || nrow(nn_order) != n || ncol(nn_order) != (n - 1))
      stop("nn_order must be an n x (n-1) integer matrix matching D.")
    if (any(nn_order < 1 | nn_order > n, na.rm = TRUE))
      stop("nn_order contains indices outside [1, n].")
    for (i in seq_len(n)) {
      nn_order[i, ] <- nn_order[i, nn_order[i, ] != i]
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
    k <- floor(pi_nn * (m - 1L))
    max(1L, min(m - 1L, as.integer(k)))
  }

  knn_prop_target <- function(i, allowed_mask, k, target_level) {
    row_ord <- NN_order[i, ]
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

      allowed <- rep(FALSE, n); allowed[pair_idx] <- TRUE

      # soft scores: p_a(i) for each i in the pair
      p_a <- purrr::map_dbl(pair_idx, ~ knn_prop_target(.x, allowed, k_use, L[a]))

      # true labels (binary)
      y <- as.integer(labels[pair_idx] == L[a])  # 1 for a, 0 for b
      pi_a <- mean(y)

      # if NA probs happen (should be rare), drop them safely
      ok <- !is.na(p_a)
      if (!any(ok)) next
      p <- p_a[ok]
      yk <- y[ok]
      pi_ok <- mean(yk)

      # --- scoring & mapping to [0,1] ---------------------------
      if (score == "ba") {
        thr <- if (decision == "posterior") 0.5 else pi_ok
        yhat_is_a <- (p >= thr)  # ties go to 'a'
        rec_a <- if (any(yk == 1)) mean(yhat_is_a[yk == 1]) else 0
        rec_b <- if (any(yk == 0)) mean(!yhat_is_a[yk == 0]) else 0
        BA    <- 0.5 * (rec_a + rec_b)
        val <- max(0, min(1, 2 * (BA - 0.5)))

      } else if (score == "brier") {
        # Brier score and normalized improvement vs prior-only predictor
        bs  <- mean((yk - p)^2)
        bs0 <- pi_ok * (1 - pi_ok)  # baseline Brier for constant predictor = prior
        val <- if (bs0 <= 0) 0 else max(0, min(1, 1 - bs / bs0))

      } else { # score == "logloss"
        # Log loss and normalized improvement vs prior-only predictor
        p2 <- pmin(pmax(p, eps), 1 - eps)
        ll  <- -mean(yk * log(p2) + (1 - yk) * log(1 - p2))
        ll0 <- -(pi_ok * log(pmin(pmax(pi_ok, eps), 1 - eps)) +
                   (1 - pi_ok) * log(pmin(pmax(1 - pi_ok, eps), 1 - eps)))
        val <- if (ll0 <= 0) 0 else max(0, min(1, 1 - ll / ll0))
      }

      Delta[a, b] <- val
      Delta[b, a] <- val
    }
  }

  diag(Delta) <- 0
  Delta
}
