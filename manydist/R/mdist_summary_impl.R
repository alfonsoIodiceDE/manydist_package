mdist_summary_impl <- function(object, ...) {
  m <- as.matrix(object$distance)
  p <- object$params
  n_row <- nrow(m); n_col <- ncol(m)
  
  cat("MDist summary\n")
  cat("  Preset :", object$preset, "\n")
  if (n_row == n_col) cat("  Observations :", n_row, "\n")
  else cat("  Train–test matrix :", n_col, "train ×", n_row, "test\n")
  
  if (!is.null(p$cont_p)) cat("  Continuous vars:", p$cont_p, "\n")
  if (!is.null(p$cat_p))  cat("  Categorical vars:", p$cat_p, "\n")
  
  cat("\nDistance matrix statistics:\n")
  cat("  Mean :", round(mean(m), 4), "\n")
  cat("  SD   :", round(sd(m), 4), "\n")
  cat("  Min  :", round(min(m), 4), "\n")
  cat("  Max  :", round(max(m), 4), "\n")
  cat("  Symmetric:", isTRUE(all.equal(m, t(m), tolerance = 1e-10)), "\n\n")
  
  has_cont <- !is.null(p$cont_p) && p$cont_p > 0
  has_cat  <- !is.null(p$cat_p)  && p$cat_p  > 0
  msg <- if (has_cont && has_cat) {
    paste0(if (isTRUE(p$commensurable)) "Commensurable" else "Non-commensurable",
           " mixed-type distance combining ",
           p$distance_cont, " (continuous) and ",
           p$distance_cat, " (categororical) components.")
  } else if (has_cont) {
    paste("Purely continuous distance using", p$distance_cont, "metric.")
  } else if (has_cat) {
    paste("Purely categorical distance using", p$distance_cat, "metric.")
  } else "Empty or untyped distance object."
  cat(msg, "\n")
  
  invisible(object)
}
