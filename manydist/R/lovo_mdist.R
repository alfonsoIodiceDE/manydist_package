# R/lovo_mdist.R

# internal helper (not exported)
.get_partition <- function(D, method, k, hclust_method = "average") {
  if (method == "pam") {
    return(cluster::pam(stats::as.dist(D), k = k, diss = TRUE)$clustering)
  }

  if (method == "hclust") {
    hc <- stats::hclust(stats::as.dist(D), method = hclust_method)
    return(stats::cutree(hc, k = k))
  }

  stop("Unknown clustering method.")
}

MDistLOVO <- R6::R6Class(
  "MDistLOVO",
  public = list(
    results         = NULL,
    base_mds        = NULL,
    dims            = NULL,
    preset          = NULL,
    params          = NULL,
    n_obs           = NULL,
    full_dist       = NULL,
    loo_dist        = NULL,
    cluster_k       = NULL,
    cluster_methods = NULL,
    hclust_method   = NULL,
    response_used   = NULL,

    initialize = function(x, response = NULL, ..., dims = 2, keep_dist = FALSE,
                          cluster_k = NULL,
                          cluster_methods = c("pam", "hclust"),
                          hclust_method = "average",
                          response_used = TRUE) {

      x <- tibble::as_tibble(x)
      dots <- rlang::list2(...)

      cluster_methods <- unique(cluster_methods)
      valid_methods <- c("pam", "hclust")

      if (!all(cluster_methods %in% valid_methods)) {
        stop("cluster_methods must be a subset of c('pam', 'hclust').")
      }

      if (!is.logical(response_used) || length(response_used) != 1L || is.na(response_used)) {
        stop("response_used must be a single TRUE/FALSE value.")
      }

      if (!is.null(cluster_k)) {
        if (!is.numeric(cluster_k) || length(cluster_k) != 1 || is.na(cluster_k)) {
          stop("cluster_k must be a single non-missing numeric value.")
        }
        cluster_k <- as.integer(cluster_k)
        if (cluster_k < 2 || cluster_k >= nrow(x)) {
          stop("cluster_k must be an integer between 2 and nrow(x) - 1.")
        }
      }

      response_col <- NULL
      if (!is.null(response)) {
        if (!is.character(response) || length(response) != 1L) {
          stop("`response` must be a single column name.", call. = FALSE)
        }
        response_col <- response
      }

      if (!is.null(response_col) && !response_col %in% names(x)) {
        stop("`response` must name a column inside `x`.", call. = FALSE)
      }

      if (isTRUE(response_used)) {
        x_work <- x
        mdist_args <- c(list(x = x_work), dots)
        if (!is.null(response_col)) {
          mdist_args$response <- response_col
        }
      } else {
        x_work <- x

        if (!is.null(response_col) && response_col %in% names(x_work)) {
          x_work <- x_work |>
            dplyr::select(-dplyr::all_of(response_col))
        }

        mdist_args <- c(list(x = x_work), dots)
      }


      md_full <- rlang::exec(mdist, !!!mdist_args)
      self$preset <- md_full$preset
      self$params <- md_full$params
      self$response_used <- response_used

      # variables eligible for LOVO are always predictors only
      x_predictors <- x_work
      if (!is.null(response_col) && response_col %in% names(x_predictors)) {
        x_predictors <- x_predictors |>
          dplyr::select(-dplyr::all_of(response_col))
      }

      var_types <- purrr::map_chr(
        .x = x_predictors,
        .f = ~ if (is.numeric(.x)) "numeric" else "categorical"
      )
      full_mat <- md_full$to_dist() |> as.matrix()
      self$n_obs <- nrow(full_mat)
      self$dims  <- dims

      self$cluster_k       <- cluster_k
      self$cluster_methods <- cluster_methods
      self$hclust_method   <- hclust_method

      if (nrow(full_mat) != ncol(full_mat)) {
        stop("cmdscale() needs a square distance. Avoid validate_x here or handle train-only MDS.")
      }

      max_dims <- nrow(full_mat) - 1
      if (dims > max_dims) {
        stop(sprintf("dims must be <= %d.", max_dims))
      }

      full_partitions <- list()

      if (!is.null(cluster_k)) {
        if ("pam" %in% cluster_methods) {
          full_partitions$pam <- .get_partition(
            D = full_mat, method = "pam", k = cluster_k,
            hclust_method = hclust_method
          )
        }
        if ("hclust" %in% cluster_methods) {
          full_partitions$hclust <- .get_partition(
            D = full_mat, method = "hclust", k = cluster_k,
            hclust_method = hclust_method
          )
        }
      }

      self$base_mds <- stats::cmdscale(full_mat, eig = TRUE, k = dims)$points[, 1:dims, drop = FALSE]

      if (keep_dist) {
        self$full_dist <- full_mat
      }

      vars <- names(x_predictors)

      loo_list <- vector("list", length(vars))
      names(loo_list) <- vars

      ari_pam <- rep(NA_real_, length(vars))
      ari_hclust <- rep(NA_real_, length(vars))

      for (i in seq_along(vars)) {
        var <- vars[i]

        x_subset <- x_work |>
          dplyr::select(-dplyr::any_of(var))

        md_loo <- rlang::exec(mdist, x = x_subset, !!!mdist_args[names(mdist_args) != "x"])
        loo_mat <- md_loo$to_dist() |> as.matrix()
        loo_list[[i]] <- loo_mat

        if (!is.null(cluster_k)) {
          if ("pam" %in% cluster_methods) {
            loo_pam <- .get_partition(
              D = loo_mat, method = "pam", k = cluster_k,
              hclust_method = hclust_method
            )
            ari_pam[i] <- mclust::adjustedRandIndex(full_partitions$pam, loo_pam)
          }

          if ("hclust" %in% cluster_methods) {
            loo_hc <- .get_partition(
              D = loo_mat, method = "hclust", k = cluster_k,
              hclust_method = hclust_method
            )
            ari_hclust[i] <- mclust::adjustedRandIndex(full_partitions$hclust, loo_hc)
          }
        }
      }

      if (keep_dist) {
        self$loo_dist <- loo_list
      }

      mad <- vapply(loo_list, function(m) mean(abs(full_mat - m)), numeric(1))

      cc <- vapply(loo_list, function(m) {
        pts <- stats::cmdscale(m, eig = TRUE, k = dims)$points[, 1:dims, drop = FALSE]
        congruence_coeff(self$base_mds, pts)
      }, numeric(1))

      ac <- sqrt(1 - cc^2)

      res <- tibble::tibble(
        variable          = vars,
        variable_type     = unname(var_types[vars]),
        mad_importance    = mad,
        cc_importance     = cc,
        mds_congruence    = cc,
        ac_importance     = ac,
        mad_normalized    = mad / sum(mad),
        relative_distance = mad / sum(mad)
      )

      if (!is.null(cluster_k)) {
        if ("pam" %in% cluster_methods) {
          res <- res |>
            dplyr::mutate(
              ari_pam = ari_pam,
              pam_importance = 1 - ari_pam
            )
        }

        if ("hclust" %in% cluster_methods) {
          res <- res |>
            dplyr::mutate(
              ari_hclust = ari_hclust,
              hclust_importance = 1 - ari_hclust
            )
        }
      }

      self$results <- res
    },

    print = function(...) {
      cat("MDistLOVO object\n")
      cat("  preset :", self$preset, "\n")
      cat("  dims   :", self$dims, "\n")
      cat("  n_obs  :", self$n_obs, "\n")
      cat("  response used :", self$response_used, "\n")

      if (!is.null(self$cluster_k)) {
        cat("  cluster_k :", self$cluster_k, "\n")
        cat("  cluster methods :", paste(self$cluster_methods, collapse = ", "), "\n")
        if ("hclust" %in% self$cluster_methods) {
          cat("  hclust linkage :", self$hclust_method, "\n")
        }
      }

      r <- self$results

      has_pam <- "pam_importance" %in% names(r) && !all(is.na(r$pam_importance))
      has_hclust <- "hclust_importance" %in% names(r) && !all(is.na(r$hclust_importance))

      if (has_pam || has_hclust) {
        cat("  clustering diagnostics :")
        if (has_pam) cat(" PAM")
        if (has_hclust) cat(" HCLUST")
        cat("\n")
      }

      ord_metric <- if ("relative_distance" %in% names(r)) "relative_distance" else "mad_importance"
      ord <- order(r[[ord_metric]], decreasing = TRUE, na.last = NA)

      top_cols <- c("variable", "variable_type")

      if ("relative_distance" %in% names(r)) {
        top_cols <- c(top_cols, "relative_distance")
      } else if ("mad_normalized" %in% names(r)) {
        top_cols <- c(top_cols, "mad_normalized")
      }

      top_cols <- c(top_cols, "mad_importance")

      if (has_pam) {
        top_cols <- c(top_cols, "pam_importance")
      }
      if (has_hclust) {
        top_cols <- c(top_cols, "hclust_importance")
      }

      top <- utils::head(r[ord, top_cols, drop = FALSE], 5)

      cat("  top vars:\n")
      if (nrow(top) == 0) {
        cat("   no finite LOVO importance values available\n")
      } else {
        print(top, row.names = FALSE)
      }

      invisible(self)
    },

    summary = function(...) {
      cat("Summary of MDistLOVO\n")
      cat("  preset :", self$preset, "\n")
      cat("  dims   :", self$dims, "\n")
      cat("  n_obs  :", self$n_obs, "\n")
      cat("  response used :", self$response_used, "\n")

      if (!is.null(self$cluster_k)) {
        cat("  cluster_k :", self$cluster_k, "\n")
        cat("  cluster methods :", paste(self$cluster_methods, collapse = ", "), "\n")
        if ("hclust" %in% self$cluster_methods) {
          cat("  hclust linkage :", self$hclust_method, "\n")
        }
      }

      cat("\n")

      r <- self$results

      .print_metric_summary <- function(x, label) {
        if (all(is.na(x))) return(invisible(NULL))
        xr <- range(x, na.rm = TRUE)
        xm <- mean(x, na.rm = TRUE)
        cat(label, ":\n", sep = "")
        cat(sprintf("  range [%.4f, %.4f], mean %.4f\n\n", xr[1], xr[2], xm))
      }

      .print_top <- function(df, col, label, n = 5, decreasing = TRUE) {
        if (!(col %in% names(df)) || all(is.na(df[[col]]))) return(invisible(NULL))
        ord <- order(df[[col]], decreasing = decreasing, na.last = NA)
        cols_to_show <- c("variable", "variable_type", col)
        cols_to_show <- cols_to_show[cols_to_show %in% names(df)]
        cat(label, ":\n", sep = "")
        print(utils::head(df[ord, cols_to_show, drop = FALSE], n), row.names = FALSE)
        cat("\n")
      }

      .print_metric_summary(r$mad_importance, "MAD importance")

      if ("relative_distance" %in% names(r)) {
        .print_metric_summary(r$relative_distance, "Relative distance")
      } else if ("mad_normalized" %in% names(r)) {
        .print_metric_summary(r$mad_normalized, "Normalized MAD importance")
      }

      if ("mds_congruence" %in% names(r)) {
        .print_metric_summary(r$mds_congruence, "MDS congruence")
      } else if ("cc_importance" %in% names(r)) {
        .print_metric_summary(r$cc_importance, "Congruence coefficient (CC)")
      }

      .print_metric_summary(r$ac_importance, "Alienation coefficient (AC)")

      if ("ari_pam" %in% names(r)) {
        .print_metric_summary(r$ari_pam, "ARI vs full PAM partition")
      }
      if ("pam_importance" %in% names(r)) {
        .print_metric_summary(r$pam_importance, "PAM importance (1 - ARI)")
      }

      if ("ari_hclust" %in% names(r)) {
        .print_metric_summary(r$ari_hclust, "ARI vs full HCLUST partition")
      }
      if ("hclust_importance" %in% names(r)) {
        .print_metric_summary(r$hclust_importance, "HCLUST importance (1 - ARI)")
      }

      if ("relative_distance" %in% names(r)) {
        .print_top(r, "relative_distance", "Top by relative distance")
      } else if ("mad_normalized" %in% names(r)) {
        .print_top(r, "mad_normalized", "Top by normalized MAD")
      }

      .print_top(r, "mad_importance", "Top by MAD")

      if ("mds_congruence" %in% names(r)) {
        .print_top(r, "mds_congruence", "Top by MDS congruence", decreasing = FALSE)
      } else if ("cc_importance" %in% names(r)) {
        .print_top(r, "cc_importance", "Top by congruence coefficient", decreasing = FALSE)
      }

      .print_top(r, "ac_importance", "Top by alienation coefficient")

      .print_top(r, "pam_importance", "Top by PAM importance")
      .print_top(r, "hclust_importance", "Top by HCLUST importance")

      invisible(self)
    },

    autoplot = function(metric = c("relative_distance", "mad_importance",
                                   "mds_congruence", "ac_importance",
                                   "cc_importance", "mad_normalized",
                                   "ari_pam", "ari_hclust",
                                   "pam_importance", "hclust_importance"),
                        reorder = FALSE,
                        top_n = NULL) {
      metric <- match.arg(metric)

      metric_label <- switch(
        metric,
        relative_distance = "Relative distance",
        mad_importance    = "MAD importance",
        mds_congruence    = "MDS congruence",
        ac_importance     = "Alienation coefficient",
        cc_importance     = "Congruence coefficient",
        mad_normalized    = "Normalized MAD importance",
        ari_pam           = "ARI vs full PAM partition",
        ari_hclust        = "ARI vs full HCLUST partition",
        pam_importance    = "PAM importance (1 - ARI)",
        hclust_importance = "HCLUST importance (1 - ARI)"
      )

      df <- self$results |>
        dplyr::mutate(method = self$preset %||% "lovo")

      if (!(metric %in% names(df))) {
        stop(sprintf("Metric '%s' not found in results.", metric))
      }

      if (all(is.na(df[[metric]]))) {
        stop(sprintf("Metric '%s' is available but contains only NA values.", metric))
      }

      smaller_is_stronger <- metric %in% c("ari_pam", "ari_hclust", "cc_importance", "mds_congruence")

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
          title = paste("LOVO:", metric_label)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

      if (metric %in% c("ari_pam", "ari_hclust",
                        "pam_importance", "hclust_importance")) {
        p <- p + ggplot2::coord_cartesian(ylim = c(0, 1))
      }

      p
    }
  )
)

# tiny factory for symmetry with mdist()
#' @export
lovo_mdist <- function(x, response = NULL, ..., dims = 2, keep_dist = FALSE,
                       cluster_k = NULL,
                       cluster_methods = c("pam", "hclust"),
                       hclust_method = "average",
                       response_used = TRUE) {

  response_name <- NULL
  response_quo <- rlang::enquo(response)

  if (!rlang::quo_is_null(response_quo)) {
    response_expr <- rlang::quo_get_expr(response_quo)

    if (rlang::is_string(response_expr)) {
      response_name <- response_expr
    } else {
      response_name <- rlang::as_name(response_expr)
    }
  }

  MDistLOVO$new(
    x = x,
    response = response_name,
    ...,
    dims = dims,
    keep_dist = keep_dist,
    cluster_k = cluster_k,
    cluster_methods = cluster_methods,
    hclust_method = hclust_method,
    response_used = response_used
  )
}
