.get_step_mdist <- function(prepped_recipe) {
  idx <- purrr::detect_index(prepped_recipe$steps, ~ inherits(.x, "step_mdist"))
  if (idx == 0) {
    stop("No `step_mdist()` found in the recipe.", call. = FALSE)
  }
  prepped_recipe$steps[[idx]]
}

.extract_dist_block <- function(x, prefix = "^dist_") {
  x <- tibble::as_tibble(x)
  dist_cols <- grep(prefix, names(x), value = TRUE)

  if (length(dist_cols) == 0L) {
    stop("No distance columns found in baked data.", call. = FALSE)
  }

  as.matrix(x[, dist_cols, drop = FALSE])
}
