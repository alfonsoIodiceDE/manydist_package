# R/lovo_mdist_compare.R

MDistLOVOCompare <- R6::R6Class(
  "MDistLOVOCompare",
  public = list(
    results = NULL,
    methods = NULL,
    dims    = NULL,
    n_obs   = NULL,

    initialize = function(results, methods, dims = 2, n_obs = NA_integer_) {
      self$results <- tibble::as_tibble(results)
      self$methods <- methods
      self$dims    <- dims
      self$n_obs   <- n_obs
    },

    print = function(...) {
      cat("MDistLOVOCompare object\n")
      cat("  methods:", paste(names(self$methods), collapse = ", "), "\n")
      cat("  dims   :", self$dims, "\n")
      cat("  n_obs  :", self$n_obs, "\n")
      cat("  rows   :", nrow(self$results), "\n\n")

      r <- self$results

      top_cols <- c("method", "variable", "mad_importance", "ac_importance")
      if ("pam_importance" %in% names(r) && !all(is.na(r$pam_importance))) {
        top_cols <- c(top_cols, "pam_importance")
      }
      if ("hclust_importance" %in% names(r) && !all(is.na(r$hclust_importance))) {
        top_cols <- c(top_cols, "hclust_importance")
      }

      top <- r |>
        dplyr::group_by(method) |>
        dplyr::slice_max(order_by = mad_importance, n = 3, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(top_cols))

      cat("Top variables by MAD within each method:\n")
      print(top, n = min(nrow(top), 12))
      invisible(self)
    },

    summary = function(...) {
      r <- self$results

      out <- r |>
        dplyr::group_by(method) |>
        dplyr::summarise(
          mad_min  = min(mad_importance, na.rm = TRUE),
          mad_max  = max(mad_importance, na.rm = TRUE),
          mad_mean = mean(mad_importance, na.rm = TRUE),
          ac_min   = min(ac_importance, na.rm = TRUE),
          ac_max   = max(ac_importance, na.rm = TRUE),
          ac_mean  = mean(ac_importance, na.rm = TRUE),
          cc_min   = min(cc_importance, na.rm = TRUE),
          cc_max   = max(cc_importance, na.rm = TRUE),
          cc_mean  = mean(cc_importance, na.rm = TRUE),
          .groups = "drop"
        )

      if ("ari_pam" %in% names(r) && !all(is.na(r$ari_pam))) {
        out <- out |>
          dplyr::left_join(
            r |>
              dplyr::group_by(method) |>
              dplyr::summarise(
                ari_pam_min  = min(ari_pam, na.rm = TRUE),
                ari_pam_max  = max(ari_pam, na.rm = TRUE),
                ari_pam_mean = mean(ari_pam, na.rm = TRUE),
                pam_imp_min  = min(pam_importance, na.rm = TRUE),
                pam_imp_max  = max(pam_importance, na.rm = TRUE),
                pam_imp_mean = mean(pam_importance, na.rm = TRUE),
                .groups = "drop"
              ),
            by = "method"
          )
      }

      if ("ari_hclust" %in% names(r) && !all(is.na(r$ari_hclust))) {
        out <- out |>
          dplyr::left_join(
            r |>
              dplyr::group_by(method) |>
              dplyr::summarise(
                ari_hclust_min  = min(ari_hclust, na.rm = TRUE),
                ari_hclust_max  = max(ari_hclust, na.rm = TRUE),
                ari_hclust_mean = mean(ari_hclust, na.rm = TRUE),
                hclust_imp_min  = min(hclust_importance, na.rm = TRUE),
                hclust_imp_max  = max(hclust_importance, na.rm = TRUE),
                hclust_imp_mean = mean(hclust_importance, na.rm = TRUE),
                .groups = "drop"
              ),
            by = "method"
          )
      }

      cat("Summary of MDistLOVOCompare\n")
      cat("  methods:", paste(names(self$methods), collapse = ", "), "\n")
      cat("  dims   :", self$dims, "\n")
      cat("  n_obs  :", self$n_obs, "\n\n")
      print(out, n = nrow(out))
      invisible(out)
    },

    autoplot = function(metric = c("mad_importance", "mad_normalized",
                                   "ac_importance", "cc_importance",
                                   "ari_pam", "ari_hclust",
                                   "pam_importance", "hclust_importance"),
                        reorder = FALSE,
                        top_n = NULL) {
      metric <- match.arg(metric)

      metric_label <- switch(
        metric,
        mad_importance    = "MAD importance",
        mad_normalized    = "Normalized MAD importance",
        ac_importance     = "Alienation coefficient importance",
        cc_importance     = "Congruence coefficient",
        ari_pam           = "ARI vs full PAM partition",
        ari_hclust        = "ARI vs full HCLUST partition",
        pam_importance    = "PAM importance (1 - ARI)",
        hclust_importance = "HCLUST importance (1 - ARI)"
      )

      df <- self$results

      if (!(metric %in% names(df))) {
        stop(sprintf("Metric '%s' not found in results.", metric))
      }

      if (all(is.na(df[[metric]]))) {
        stop(sprintf("Metric '%s' is available but contains only NA values.", metric))
      }

      smaller_is_stronger <- metric %in% c("ari_pam", "ari_hclust", "cc_importance")

      if (!is.null(top_n)) {
        top_vars <- df |>
          dplyr::group_by(variable) |>
          dplyr::summarise(
            avg_metric = mean(.data[[metric]], na.rm = TRUE),
            .groups = "drop"
          )

        if (smaller_is_stronger) {
          top_vars <- top_vars |>
            dplyr::slice_min(order_by = avg_metric, n = top_n, with_ties = FALSE)
        } else {
          top_vars <- top_vars |>
            dplyr::slice_max(order_by = avg_metric, n = top_n, with_ties = FALSE)
        }

        top_vars <- top_vars |>
          dplyr::pull(variable)

        df <- df |>
          dplyr::filter(variable %in% top_vars)
      }

      if (reorder) {
        ord_df <- df |>
          dplyr::group_by(variable) |>
          dplyr::summarise(
            avg_metric = mean(.data[[metric]], na.rm = TRUE),
            .groups = "drop"
          )

        if (smaller_is_stronger) {
          ord <- ord_df |>
            dplyr::arrange(avg_metric) |>
            dplyr::pull(variable)
        } else {
          ord <- ord_df |>
            dplyr::arrange(dplyr::desc(avg_metric)) |>
            dplyr::pull(variable)
        }
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
            alpha = 0.08
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
              alpha = 0.08
            )
        }
      }

      p <- p +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(
          x = "left-out variable",
          y = metric_label,
          color = "method",
          title = paste("LOVO comparison:", metric_label)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

      if (metric %in% c("mad_normalized", "ac_importance",
                        "ari_pam", "ari_hclust",
                        "pam_importance", "hclust_importance")) {
        p <- p + ggplot2::coord_cartesian(ylim = c(0, 1))
      }

      p
    }
  )
)


# -----------------------------------------------------------
# User-facing constructor
# -----------------------------------------------------------

#' Compare LOVO diagnostics across multiple distance specifications
#'
#' This function compares **leave-one-variable-out (LOVO)** diagnostics
#' across multiple distance definitions supported by \code{manydist}.
#'
#' For each distance specification, the function:
#'
#' \enumerate{
#' \item Computes the full mixed-type distance using \code{mdist()}.
#' \item Recomputes the distance repeatedly leaving out one variable at a time.
#' \item Measures the impact of each variable using metrics such as
#' mean absolute deviation (MAD), congruence-based diagnostics, and,
#' when requested, clustering-based agreement measures.
#' }
#'
#' The results are combined across methods and returned as an
#' \code{MDistLOVOCompare} object, which supports
#' \code{print()}, \code{summary()}, and \code{ggplot2::autoplot()}.
#'
#' @param x A data frame or tibble containing the predictors.
#'
#' @param methods A **named list** describing the distance specifications
#' to compare. Each element must be a list of arguments passed to
#' \code{\link{lovo_mdist}}.
#'
#' For example:
#'
#' \preformatted{
#' methods = list(
#'   gower = list(preset = "gower"),
#'   u_dep = list(preset = "unbiased_dependent"),
#'   custom = list(
#'     distance_cont = "manhattan",
#'     distance_cat  = "matching",
#'     commensurable = TRUE
#'   )
#' )
#' }
#'
#' @param dims Number of dimensions used for the MDS configuration when
#' computing congruence-based diagnostics.
#'
#' @param keep_dist Logical; if \code{TRUE}, distance matrices from the
#' LOVO computations are retained. This increases memory usage.
#'
#' @param .progress Logical; if \code{TRUE}, progress messages are printed
#' while computing LOVO diagnostics for each method.
#'
#' @param ... Additional arguments passed to \code{\link{lovo_mdist}}
#' unless overridden in the method-specific argument list.
#'
#' These may include optional clustering diagnostics, for example
#' \code{cluster_k}, \code{cluster_methods}, and \code{hclust_method}.
#'
#' @return An object of class \code{MDistLOVOCompare} containing:
#'
#' \describe{
#' \item{results}{A tibble with one row per method-variable combination.}
#' \item{methods}{The list of distance specifications used.}
#' \item{dims}{Number of MDS dimensions used.}
#' \item{n_obs}{Number of observations in the dataset.}
#' }
#'
#' @seealso
#' \code{\link{lovo_mdist}}, \code{\link{mdist}}
#'
#' @examples
#' \dontrun{
#' library(manydist)
#' library(palmerpenguins)
#'
#' data <- penguins |>
#'   dplyr::select(-species) |>
#'   tidyr::drop_na()
#'
#' cmp <- compare_lovo_mdist(
#'   x = data,
#'   methods = list(
#'     gower = list(preset = "gower"),
#'     u_dep = list(preset = "unbiased_dependent")
#'   ),
#'   cluster_k = 3
#' )
#'
#' summary(cmp)
#' autoplot(cmp, metric = "mad_importance")
#' autoplot(cmp, metric = "pam_importance")
#' }
#'
#' @export
compare_lovo_mdist <- function(x,
                               methods,
                               dims = 2,
                               keep_dist = FALSE,
                               .progress = FALSE,
                               ...) {
  if (!is.list(methods) || is.null(names(methods)) || any(names(methods) == "")) {
    stop("`methods` must be a named list.")
  }

  x <- tibble::as_tibble(x)

  var_types <- vapply(x, function(col) {
    if (is.numeric(col)) "numeric"
    else "categorical"
  }, character(1))

  res <- purrr::imap_dfr(methods, function(method_args, method_name) {
    if (.progress) {
      message("Running LOVO for method: ", method_name)
    }

    if (!is.list(method_args)) {
      stop("Each element of `methods` must be a list of arguments.")
    }

    args <- c(
      list(x = x, dims = dims, keep_dist = keep_dist),
      list(...),
      method_args
    )

    obj <- do.call(lovo_mdist, args)

    obj$results |>
      dplyr::mutate(
        method = method_name,
        variable_type = unname(var_types[variable]),
        .before = 1
      )
  }) |>
    dplyr::mutate(method = factor(method, levels = names(methods)))

  MDistLOVOCompare$new(
    results = res,
    methods = methods,
    dims = dims,
    n_obs = nrow(x)
  )
}
