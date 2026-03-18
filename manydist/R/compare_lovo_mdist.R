# R/lovo_mdist_compare.R

MDistLOVOCompare <- R6::R6Class(
  "MDistLOVOCompare",
  public = list(
    results = NULL,   # tibble: method, variable, variable_type, mad_importance, cc_importance, ac_importance, mad_normalized
    methods = NULL,   # named list of method specs
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

      top <- self$results |>
        dplyr::group_by(method) |>
        dplyr::slice_max(order_by = mad_importance, n = 3, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(method, variable, mad_importance, ac_importance)

      cat("Top variables by MAD within each method:\n")
      print(top, n = min(nrow(top), 12))
      invisible(self)
    },

    summary = function(...) {
      out <- self$results |>
        dplyr::group_by(method) |>
        dplyr::summarise(
          mad_min  = min(mad_importance, na.rm = TRUE),
          mad_max  = max(mad_importance, na.rm = TRUE),
          mad_mean = mean(mad_importance, na.rm = TRUE),
          ac_min   = min(ac_importance, na.rm = TRUE),
          ac_max   = max(ac_importance, na.rm = TRUE),
          ac_mean  = mean(ac_importance, na.rm = TRUE),
          .groups = "drop"
        )

      cat("Summary of MDistLOVOCompare\n")
      cat("  methods:", paste(names(self$methods), collapse = ", "), "\n")
      cat("  dims   :", self$dims, "\n")
      cat("  n_obs  :", self$n_obs, "\n\n")
      print(out, n = nrow(out))
      invisible(out)
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

      df <- self$results

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
          title = paste("LOVO comparison:", metric_label)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
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
#' mean absolute deviation (MAD) and the alienation coefficient.
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
#'   )
#' )
#'
#' summary(cmp)
#' autoplot(cmp, metric = "mad_importance")
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
