#' Internal: summary method implementation for MDist R6 objects
mdist_summary_impl <- function(object, ...) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  d <- object$distance
  if (is.null(d)) {
    cat("MDist summary\n  Preset :", object$preset, "\n  <no distance stored>\n")
    return(invisible(object))
  }
  
  is_dist <- inherits(d, "dist")
  
  # Dimensions without full materialization
  if (is_dist) {
    n <- attr(d, "Size")
    n_row <- n_col <- n
    is_square <- TRUE
  } else {
    dm <- dim(d)
    n_row <- dm[1]; n_col <- dm[2]
    is_square <- (n_row == n_col)
  }
  
  # Helper: fetch a tiny block (≤10×10) for cheap diagnostics/preview
  get_block <- function(x, k = 10) {
    if (inherits(x, "dist")) {
      idx <- seq_len(min(k, attr(x, "Size")))
      as.matrix(x)[idx, idx, drop = FALSE]
    } else {
      r <- seq_len(min(k, nrow(x)))
      c <- seq_len(min(k, ncol(x)))
      as.matrix(x[r, c, drop = FALSE])
    }
  }
  block <- get_block(d)
  
  # Stats (exclude diagonal for square matrices). For dist, use the vector form directly.
  if (is_dist) {
    dv <- as.numeric(d)
    rng <- range(dv, na.rm = TRUE)
    m_mean <- mean(dv, na.rm = TRUE)
    m_sd   <- stats::sd(dv, na.rm = TRUE)
    na_cnt <- sum(is.na(dv))
    inf_cnt <- sum(is.infinite(dv))
  } else {
    # For rectangular or huge square matrices, compute stats on the small block as a proxy
    # (full scan would be too costly). This keeps the function fast and memory-safe.
    b <- block
    if (is_square) diag(b) <- NA_real_
    rng <- range(b, na.rm = TRUE)
    m_mean <- mean(b, na.rm = TRUE)
    m_sd   <- stats::sd(b, na.rm = TRUE)
    na_cnt <- sum(is.na(b))
    inf_cnt <- sum(is.infinite(b))
  }
  
  # Cheap symmetry & diagonal checks (on block only when applicable)
  diag_zero <- if (is_square) all(diag(block) == 0) else NA
  approx_sym <- if (is_square) {
    max_abs_diff <- max(abs(block - t(block)), na.rm = TRUE)
    max_abs_diff < 1e-12
  } else NA
  
  # Optional: quick-and-dirty triangle inequality check on the block
  tri_viol_rate <- if (is_square && nrow(block) >= 3) {
    k <- nrow(block)
    I <- sample.int(k, min(50, k), replace = TRUE)
    J <- sample.int(k, min(50, k), replace = TRUE)
    L <- sample.int(k, min(50, k), replace = TRUE)
    v <- 0L; tot <- 0L
    for (t in seq_along(I)) {
      i <- I[t]; j <- J[t]; l <- L[t]
      if (i != j && j != l && i != l) {
        tot <- tot + 1L
        if (!(block[i, j] <= block[i, l] + block[l, j])) v <- v + 1L
      }
    }
    if (tot > 0) v / tot else NA_real_
  } else NA_real_
  
  p <- object$params %||% list()
  has_cont <- !is.null(p$cont_p) && isTRUE(p$cont_p > 0)
  has_cat  <- !is.null(p$cat_p)  && isTRUE(p$cat_p  > 0)
  
  cat("MDist summary\n")
  cat("  Preset :", object$preset, "\n")
  if (is_square) cat("  Observations :", n_row, "\n")
  else           cat("  Train–test matrix :", n_col, "train ×", n_row, "test\n")
  
  if (!is.null(p$cont_p)) cat("  Continuous vars :", p$cont_p, "\n")
  if (!is.null(p$cat_p))  cat("  Categorical vars:", p$cat_p, "\n")
  
  cat("\nDistance matrix statistics",
      if (!is_dist && is_square) " (on preview block)" else if (!is_dist) " (on preview block)" else "",
      ":\n", sep = "")
  cat("  Mean :", round(m_mean, 4), "\n")
  cat("  SD   :", round(m_sd,   4), "\n")
  cat("  Min  :", round(rng[1], 4), "\n")
  cat("  Max  :", round(rng[2], 4), "\n")
  if (is_square) {
    cat("  Diagonal all zeros (block):", if (isTRUE(diag_zero)) "YES" else "NO", "\n")
    cat("  Symmetric (block check)   :", if (isTRUE(approx_sym)) "YES" else "NO", "\n")
    if (!is.na(tri_viol_rate))
      cat("  Triangle viol. rate (block):", sprintf("%.3f", tri_viol_rate), "\n")
  }
  cat("  NA count",
      if (is_dist) "" else " (block)",
      ":", na_cnt, "\n", sep = "")
  cat("  Inf count",
      if (is_dist) "" else " (block)",
      ":", inf_cnt, "\n\n", sep = "")
  
  msg <- if (has_cont && has_cat) {
    paste0(if (isTRUE(p$commensurable)) "Commensurable" else "Non-commensurable",
           " mixed-type distance combining ",
           p$distance_cont %||% "<unspecified>", " (continuous) and ",
           p$distance_cat  %||% "<unspecified>", " (categorical) components.")
  } else if (has_cont) {
    paste("Purely continuous distance using", p$distance_cont %||% "<unspecified>", "metric.")
  } else if (has_cat) {
    paste("Purely categorical distance using", p$distance_cat %||% "<unspecified>", "metric.")
  } else "Empty or untyped distance object."
  cat(msg, "\n")
  
  invisible(object)
}
