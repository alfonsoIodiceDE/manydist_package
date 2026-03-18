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

      var_types <- vapply(x, function(col) {
        if (is.numeric(col)) "numeric"
        else "categorical"
      }, character(1))

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
        variable_type  = unname(var_types[vars]),
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

    autoplot = function(metric = c("mad_importance", "mad_normalized",
                                   "ac_importance", "cc_importance"),
                        reorder = FALSE,
                        top_n = NULL) {
      metric <- match.arg(metric)

      metric_label <- switch(
        metric,
        mad_importance = "MAD importance",
        mad_normalized = "Normalized MAD importance",
        ac_importance  = "Alienation coefficient importance",
        cc_importance  = "Congruence coefficient"
      )

      df <- self$results |>
        dplyr::mutate(
          method = self$preset %||% "lovo"
        )

      if (!is.null(top_n)) {
        top_vars <- df |>
          dplyr::group_by(variable) |>
          dplyr::summarise(
            avg_metric = mean(.data[[metric]], na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::slice_max(order_by = avg_metric, n = top_n, with_ties = FALSE) |>
          dplyr::pull(variable)

        df <- df |>
          dplyr::filter(variable %in% top_vars)
      }

      if (reorder) {
        ord <- df |>
          dplyr::group_by(variable) |>
          dplyr::summarise(
            avg_metric = mean(.data[[metric]], na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::arrange(dplyr::desc(avg_metric)) |>
          dplyr::pull(variable)
      } else if ("variable_type" %in% names(df)) {
        cat_vars <- df |>
          dplyr::filter(variable_type == "categorical") |>
          dplyr::pull(variable) |>
          unique()

        num_vars <- df |>
          dplyr::filter(variable_type == "numeric") |>
          dplyr::pull(variable) |>
          unique()

        ord <- c(cat_vars, num_vars)
      } else {
        ord <- df |>
          dplyr::pull(variable) |>
          unique()
      }

      df <- df |>
        dplyr::mutate(variable = factor(variable, levels = ord))

      rect_data <- NULL

      if ("variable_type" %in% names(df)) {
        cats <- df |>
          dplyr::filter(variable_type == "categorical") |>
          dplyr::pull(variable) |>
          unique()

        n_cat <- length(cats)
        n_var <- length(levels(df$variable))

        if (n_cat > 0 && n_cat < n_var) {
          rect_data <- tibble::tibble(
            xmin_cat = 0.5,
            xmax_cat = n_cat + 0.5,
            xmin_num = n_cat + 0.5,
            xmax_num = n_var + 0.5
          )
        } else if (n_cat > 0) {
          rect_data <- tibble::tibble(
            xmin_cat = 0.5,
            xmax_cat = n_cat + 0.5,
            xmin_num = NA_real_,
            xmax_num = NA_real_
          )
        }
      }

      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = variable,
          y = .data[[metric]],
          color = method,
          group = method
        )
      )

      if (!is.null(rect_data)) {
        p <- p +
          ggplot2::geom_rect(
            data = rect_data,
            ggplot2::aes(
              xmin = xmin_cat,
              xmax = xmax_cat,
              ymin = -Inf,
              ymax = Inf
            ),
            inherit.aes = FALSE,
            fill = "dodgerblue",
            alpha = .08
          )

        if (!is.na(rect_data$xmin_num)) {
          p <- p +
            ggplot2::geom_rect(
              data = rect_data,
              ggplot2::aes(
                xmin = xmin_num,
                xmax = xmax_num,
                ymin = -Inf,
                ymax = Inf
              ),
              inherit.aes = FALSE,
              fill = "indianred",
              alpha = .08
            )
        }
      }

      p +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(
          x = "left-out variable",
          y = metric_label,
          color = "method",
          title = paste("LOVO:", metric_label)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )+
        ggplot2::coord_cartesian(ylim = c(0, 1))
    }
  )
)

# tiny factory for symmetry with mdist()
#' @export
lovo_mdist <- function(x, ..., dims = 2, keep_dist = FALSE) {
  MDistLOVO$new(x = x, ..., dims = dims, keep_dist = keep_dist)
}
