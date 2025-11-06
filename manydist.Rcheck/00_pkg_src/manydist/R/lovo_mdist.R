# R/lovo_mdist.R
MDistLOVO <- R6::R6Class(
  "MDistLOVO",
  public = list(
    results   = NULL,   # tibble with variable, mad, cc, ac, mad_normalized
    base_mds  = NULL,   # n x dims matrix from full distance
    dims      = NULL,
    preset    = NULL,   # inherited from md_full
    params    = NULL,   # inherited from md_full
    n_obs     = NULL,
    full_dist = NULL,   # optional (if keep_dist = TRUE)
    loo_dist  = NULL,   # optional list of matrices
    
    initialize = function(x, ..., dims = 2, keep_dist = FALSE) {
      x <- tibble::as_tibble(x)
      
      # full distance via MDist R6 → dist → matrix
      md_full <- mdist(x = x, ...)
      self$preset <- md_full$preset
      self$params <- md_full$params
      
      full_mat <- md_full$to_dist() |> as.matrix()
      self$n_obs <- nrow(full_mat)
      self$dims  <- dims
      
      if (nrow(full_mat) != ncol(full_mat)) {
        stop("cmdscale() needs a square distance. Avoid validate_x here or handle train-only MDS.")
      }
      
      self$base_mds <- cmdscale(full_mat, eig = TRUE, k = dims)$points[, 1:dims, drop = FALSE]
      if (keep_dist) self$full_dist <- full_mat
      
      vars <- names(x)
      loo_list <- vector("list", length(vars)); names(loo_list) <- vars
      
      for (i in seq_along(vars)) {
        var <- vars[i]
        x_subset <- dplyr::select(x, -dplyr::any_of(var))
        md_loo <- mdist(x = x_subset, ...)
        loo_list[[i]] <- md_loo$to_dist() |> as.matrix()
      }
      if (keep_dist) self$loo_dist <- loo_list
      
      # metrics
      mad <- vapply(loo_list, function(m) mean(abs(full_mat - m)), numeric(1))
      cc  <- vapply(loo_list, function(m) {
        pts <- cmdscale(m, eig = TRUE, k = dims)$points[, 1:dims, drop = FALSE]
        congruence_coeff(self$base_mds, pts)
      }, numeric(1))
      ac  <- sqrt(1 - cc^2)
      
      self$results <- tibble::tibble(
        variable       = vars,
        mad_importance = mad,
        cc_importance  = cc,
        ac_importance  = ac,
        mad_normalized = mad / sum(mad)
      )
    },
    
    print = function(...) {
      cat("MDistLOVO object\n")
      cat("  preset :", self$preset, "\n")
      cat("  dims   :", self$dims, "\n")
      cat("  n_obs  :", self$n_obs, "\n")
      ord <- order(self$results$mad_importance, decreasing = TRUE)
      top <- utils::head(self$results[ord, c("variable","mad_importance","ac_importance")], 5)
      cat("  top vars (MAD / AC):\n")
      print(top, row.names = FALSE)
      invisible(self)
    },
    
    summary = function(...) {
      cat("Summary of MDistLOVO\n")
      cat("  preset :", self$preset, "\n")
      cat("  dims   :", self$dims, "\n")
      cat("  n_obs  :", self$n_obs, "\n\n")
      
      r <- self$results
      mad_rng <- range(r$mad_importance, na.rm = TRUE)
      ac_rng  <- range(r$ac_importance,  na.rm = TRUE)
      
      cat("MAD:\n")
      cat(sprintf("  range [%.4f, %.4f], mean %.4f\n",
                  mad_rng[1], mad_rng[2], mean(r$mad_importance, na.rm = TRUE)))
      cat("Alienation Coefficient (AC):\n")
      cat(sprintf("  range [%.4f, %.4f], mean %.4f\n\n",
                  ac_rng[1], ac_rng[2], mean(r$ac_importance, na.rm = TRUE)))
      
      ord_mad <- order(r$mad_importance, decreasing = TRUE)
      ord_ac  <- order(r$ac_importance,  decreasing = TRUE)
      cat("Top by MAD:\n")
      print(utils::head(r[ord_mad, c("variable","mad_importance")], 5), row.names = FALSE)
      cat("\nTop by AC:\n")
      print(utils::head(r[ord_ac, c("variable","ac_importance")], 5), row.names = FALSE)
      
      invisible(self)
    },
    
    autoplot = function(metric = c("mad_normalized","ac_importance","cc_importance"), n = 15) {
      metric <- match.arg(metric)
      df <- self$results |>
        dplyr::arrange(dplyr::desc(.data[[metric]])) |>
        dplyr::slice(1:n) |>
        dplyr::mutate(variable = factor(variable, levels = rev(variable)))
      ggplot2::ggplot(df, ggplot2::aes(variable, .data[[metric]])) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = metric, title = paste("LOVO -", metric))
    }
  )
)

# tiny factory for symmetry with mdist()
#' @export
lovo_mdist <- function(x, ..., dims = 2, keep_dist = FALSE) {
  MDistLOVO$new(x = x, ..., dims = dims, keep_dist = keep_dist)
}
