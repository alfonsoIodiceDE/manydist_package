.benchmark_method_labels <- function(specs) {
  labels <- rep(NA_character_, nrow(specs))

  if ("label" %in% names(specs)) {
    supplied <- as.character(specs$label)
    use_supplied <- !is.na(supplied) & nzchar(supplied)
    labels[use_supplied] <- supplied[use_supplied]
  }

  needs_label <- is.na(labels)
  is_preset <- specs$spec_type == "preset"

  labels[needs_label & is_preset] <- specs$preset[needs_label & is_preset]

  component_rows <- needs_label & !is_preset
  labels[component_rows] <- paste0(
    specs$method_cat[component_rows],
    " + ",
    specs$method_num[component_rows],
    ifelse(
      specs$commensurable[component_rows],
      " (commensurable)",
      " (not commensurable)"
    )
  )

  make.unique(labels, sep = " #")
}

.empty_benchmark_comparisons <- function(cluster_methods = character()) {
  out <- tibble::tibble(
    method_1_id = integer(),
    method_2_id = integer(),
    method_1 = character(),
    method_2 = character(),
    mad = double(),
    relative_distance = double(),
    mds_congruence = double(),
    alienation = double()
  )

  for (method in cluster_methods) {
    out[[paste0("ari_", method)]] <- double()
  }

  out
}

.benchmark_mds <- function(D, dims) {
  tryCatch(
    stats::cmdscale(D, eig = TRUE, k = dims)$points[
      , seq_len(dims), drop = FALSE
    ],
    error = function(e) e
  )
}

.benchmark_partition <- function(D, method, cluster_k, hclust_method,
                                 spectral_sigma, spectral_nstart) {
  tryCatch(
    .get_partition(
      D = D,
      method = method,
      k = cluster_k,
      hclust_method = hclust_method,
      spectral_sigma = spectral_sigma,
      spectral_nstart = spectral_nstart
    ),
    error = function(e) e
  )
}

.benchmark_pairwise_comparisons <- function(
    results,
    successful_rows,
    labels,
    dims,
    cluster_k,
    cluster_methods,
    hclust_method,
    spectral_sigma,
    spectral_nstart
) {
  if (length(successful_rows) < 2L) {
    return(.empty_benchmark_comparisons(
      if (is.null(cluster_k)) character() else cluster_methods
    ))
  }

  matrices <- lapply(
    results[successful_rows],
    function(x) x$to_dist() |> as.matrix()
  )

  matrix_dims <- vapply(
    matrices,
    function(D) paste(dim(D), collapse = "x"),
    character(1)
  )

  if (length(unique(matrix_dims)) != 1L ||
      any(vapply(matrices, nrow, integer(1)) !=
          vapply(matrices, ncol, integer(1)))) {
    stop(
      "Successful benchmark specifications must produce square distance ",
      "matrices with identical dimensions.",
      call. = FALSE
    )
  }

  distance_vectors <- lapply(
    matrices,
    function(D) as.numeric(stats::as.dist(D))
  )
  mds_configurations <- lapply(matrices, .benchmark_mds, dims = dims)

  pairs <- utils::combn(seq_along(successful_rows), 2L)

  mad <- vapply(
    seq_len(ncol(pairs)),
    function(i) {
      first <- distance_vectors[[pairs[1L, i]]]
      second <- distance_vectors[[pairs[2L, i]]]
      valid <- is.finite(first) & is.finite(second)

      if (!any(valid)) return(NA_real_)

      mean(abs(first[valid] - second[valid]))
    },
    numeric(1)
  )

  relative_distance <- vapply(
    seq_len(ncol(pairs)),
    function(i) {
      first <- distance_vectors[[pairs[1L, i]]]
      second <- distance_vectors[[pairs[2L, i]]]
      valid <- is.finite(first) & is.finite(second)

      if (!any(valid)) return(NA_real_)

      denominator <- (
        mean(abs(first[valid])) + mean(abs(second[valid]))
      ) / 2

      if (denominator == 0) {
        return(if (mad[i] == 0) 0 else NA_real_)
      }

      mad[i] / denominator
    },
    numeric(1)
  )

  mds_congruence <- vapply(
    seq_len(ncol(pairs)),
    function(i) {
      first <- mds_configurations[[pairs[1L, i]]]
      second <- mds_configurations[[pairs[2L, i]]]

      if (inherits(first, "error") || inherits(second, "error")) {
        return(NA_real_)
      }

      value <- congruence_coeff(first, second)
      if (is.finite(value)) value else NA_real_
    },
    numeric(1)
  )

  alienation <- vapply(
    mds_congruence,
    function(value) {
      if (is.na(value)) return(NA_real_)
      sqrt(max(0, 1 - value^2))
    },
    numeric(1)
  )

  comparisons <- tibble::tibble(
    method_1_id = successful_rows[pairs[1L, ]],
    method_2_id = successful_rows[pairs[2L, ]],
    method_1 = labels[successful_rows[pairs[1L, ]]],
    method_2 = labels[successful_rows[pairs[2L, ]]],
    mad = mad,
    relative_distance = relative_distance,
    mds_congruence = mds_congruence,
    alienation = alienation
  )

  if (is.null(cluster_k)) {
    return(comparisons)
  }

  partitions <- stats::setNames(
    lapply(
      cluster_methods,
      function(method) {
        lapply(
          matrices,
          .benchmark_partition,
          method = method,
          cluster_k = cluster_k,
          hclust_method = hclust_method,
          spectral_sigma = spectral_sigma,
          spectral_nstart = spectral_nstart
        )
      }
    ),
    cluster_methods
  )

  for (method in cluster_methods) {
    method_partitions <- partitions[[method]]

    comparisons[[paste0("ari_", method)]] <- vapply(
      seq_len(ncol(pairs)),
      function(i) {
        first <- method_partitions[[pairs[1L, i]]]
        second <- method_partitions[[pairs[2L, i]]]

        if (inherits(first, "error") || inherits(second, "error")) {
          return(NA_real_)
        }

        aricode::ARI(first, second)
      },
      numeric(1)
    )
  }

  comparisons
}

#' Benchmark and compare multiple `mdist()` specifications
#'
#' Applies [mdist()] repeatedly over a tibble of distance-method
#' specifications, typically generated with [all_dist_method_specs()]. It then
#' compares every pair of successful distance specifications using distance,
#' configuration, and optional clustering diagnostics.
#'
#' Each row of `specs` is interpreted as one valid `mdist()` configuration.
#' Preset-based and custom component-based specifications are both supported.
#' Failed specifications are caught and returned in the output rather than
#' stopping the full benchmark.
#'
#' @param x A data frame or tibble of predictors, optionally including the
#'   response column.
#' @param response Optional response column inside `x`, supplied either
#'   unquoted or as a character string.
#' @param specs A tibble of method specifications. By default, this is generated
#'   with [all_dist_method_specs()]. It must contain the columns `spec_type`,
#'   `preset`, `method_cat`, `method_num`, and `commensurable`. An optional
#'   `label` column supplies display labels for comparisons and plots.
#' @param dims Integer. Number of dimensions used by classical multidimensional
#'   scaling when computing congruence and alienation coefficients.
#' @param cluster_k Optional integer. Number of clusters used for pairwise
#'   adjusted Rand indices. If `NULL`, no clustering is performed.
#' @param cluster_methods Character vector specifying the clustering methods
#'   used when `cluster_k` is supplied. Possible values are `"pam"`,
#'   `"hclust"`, and `"spectral"`.
#' @param hclust_method Character string specifying the linkage method passed
#'   to [stats::hclust()] when `"hclust"` is requested.
#' @param spectral_sigma Optional numeric value for the Gaussian affinity
#'   bandwidth used by spectral clustering. If `NULL`, the default used by
#'   [spectral_dist()] is applied.
#' @param spectral_nstart Integer. Number of random starts used by the k-means
#'   step in spectral clustering.
#'
#' @return An object of class `"MDistBenchmark"`, which is also a tibble. It
#'   contains the supplied specifications together with:
#' \describe{
#'   \item{result}{The corresponding output of [mdist()], or an error object if
#'   the specification failed.}
#'   \item{ok}{Logical indicator; `TRUE` if the run completed successfully,
#'   `FALSE` otherwise.}
#'   \item{error}{Error message for failed runs, `NA` otherwise.}
#' }
#'
#' Use [benchmark_comparisons()] to obtain the pairwise diagnostics and
#' [ggplot2::autoplot()] to draw an annotated triangular heatmap.
#'
#' @details
#' Preset specifications use the `preset` column and ignore `method_cat`,
#' `method_num`, and `commensurable`. Component specifications are evaluated as
#' `preset = "custom"` and use `method_cat`, `method_num`, and
#' `commensurable`.
#'
#' Pairwise mean absolute difference (`mad`) is computed from the lower
#' triangle of each dissimilarity matrix. The symmetric relative distance is
#' defined as
#' \deqn{
#'   \frac{2\,\mathrm{mean}(|d_a-d_b|)}
#'        {\mathrm{mean}(|d_a|)+\mathrm{mean}(|d_b|)}.
#' }
#' Classical multidimensional scaling configurations are compared using
#' [congruence_coeff()]. The corresponding alienation coefficient is
#' \eqn{\sqrt{1-c^2}}.
#'
#' When `cluster_k` is supplied, each requested clustering method is applied
#' once to every successful distance specification. Partitions produced by the
#' same clustering method are then compared pairwise using the adjusted Rand
#' index. When `cluster_k = NULL`, no clustering is performed and no ARI
#' columns are included in the comparisons.
#'
#' @examples
#' if (requireNamespace("palmerpenguins", quietly = TRUE)) {
#'   data("penguins", package = "palmerpenguins")
#'
#'   penguins_small <- palmerpenguins::penguins |>
#'     dplyr::select(
#'       species, bill_length_mm, bill_depth_mm, flipper_length_mm,
#'       body_mass_g, island, sex
#'     ) |>
#'     tidyr::drop_na()
#'
#'   specs <- all_dist_method_specs(
#'     mode = "presets_only",
#'     preset = c("gower", "u_indep", "u_dep")
#'   )
#'
#'   res <- benchmark_mdist(
#'     penguins_small,
#'     response = species,
#'     specs = specs
#'   )
#'
#'   res |>
#'     dplyr::select(spec_type, preset, ok, error)
#'
#'   benchmark_comparisons(res)
#'   ggplot2::autoplot(res, metric = "relative_distance")
#' }
#'
#' @export
benchmark_mdist <- function(
    x,
    response = NULL,
    specs = all_dist_method_specs(),
    dims = 2,
    cluster_k = NULL,
    cluster_methods = c("pam", "hclust", "spectral"),
    hclust_method = "average",
    spectral_sigma = NULL,
    spectral_nstart = 50
) {
  x <- tibble::as_tibble(x)
  specs <- tibble::as_tibble(specs)

  needed <- c(
    "spec_type", "preset",
    "method_cat", "method_num", "commensurable"
  )

  missing_cols <- setdiff(needed, names(specs))
  if (length(missing_cols) > 0) {
    stop(
      "`specs` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(dims) || length(dims) != 1L || is.na(dims) ||
      dims != as.integer(dims) || dims < 1L ||
      dims >= nrow(x)) {
    stop(
      "`dims` must be a single integer between 1 and nrow(x) - 1.",
      call. = FALSE
    )
  }
  dims <- as.integer(dims)

  cluster_methods <- unique(cluster_methods)
  valid_methods <- c("pam", "hclust", "spectral")

  if (!is.character(cluster_methods) ||
      !all(cluster_methods %in% valid_methods)) {
    stop(
      "`cluster_methods` must be a subset of ",
      "c('pam', 'hclust', 'spectral').",
      call. = FALSE
    )
  }

  if (!is.null(cluster_k)) {
    if (!is.numeric(cluster_k) || length(cluster_k) != 1L ||
        is.na(cluster_k) || cluster_k != as.integer(cluster_k)) {
      stop("`cluster_k` must be a single integer or NULL.", call. = FALSE)
    }
    cluster_k <- as.integer(cluster_k)

    if (cluster_k < 2L || cluster_k >= nrow(x)) {
      stop(
        "`cluster_k` must be an integer between 2 and nrow(x) - 1.",
        call. = FALSE
      )
    }

    if (length(cluster_methods) == 0L) {
      stop(
        "`cluster_methods` must contain at least one method when ",
        "`cluster_k` is supplied.",
        call. = FALSE
      )
    }
  }

  response_name <- NULL
  response_quo <- rlang::enquo(response)

  if (!rlang::quo_is_null(response_quo)) {
    response_expr <- rlang::quo_get_expr(response_quo)

    if (rlang::is_string(response_expr)) {
      response_name <- response_expr
    } else {
      response_name <- rlang::as_name(response_expr)
    }

    if (!response_name %in% names(x)) {
      stop("`response` must name a column inside `x`.", call. = FALSE)
    }
  }

  results <- purrr::pmap(
    specs[needed],
    function(spec_type, preset, method_cat, method_num, commensurable) {
      tryCatch(
        {
          if (identical(spec_type, "preset")) {
            if (is.null(response_name)) {
              mdist(
                x = x,
                preset = preset
              )
            } else {
              mdist(
                x = x,
                response = response_name,
                preset = preset
              )
            }
          } else {
            if (is.null(response_name)) {
              mdist(
                x = x,
                preset = "custom",
                method_cat = method_cat,
                method_num = method_num,
                commensurable = commensurable
              )
            } else {
              mdist(
                x = x,
                response = response_name,
                preset = "custom",
                method_cat = method_cat,
                method_num = method_num,
                commensurable = commensurable
              )
            }
          }
        },
        error = function(e) e
      )
    }
  )

  out <- specs |>
    dplyr::mutate(
      result = results,
      ok = !purrr::map_lgl(.data$result, inherits, what = "error"),
      error = purrr::map_chr(
        .data$result,
        ~ if (inherits(.x, "error")) conditionMessage(.x) else NA_character_
      )
    )

  labels <- .benchmark_method_labels(specs)
  successful_rows <- which(out$ok)

  comparisons <- .benchmark_pairwise_comparisons(
    results = results,
    successful_rows = successful_rows,
    labels = labels,
    dims = dims,
    cluster_k = cluster_k,
    cluster_methods = cluster_methods,
    hclust_method = hclust_method,
    spectral_sigma = spectral_sigma,
    spectral_nstart = spectral_nstart
  )

  attr(out, "comparisons") <- comparisons
  attr(out, "method_labels") <- labels
  attr(out, "benchmark_parameters") <- list(
    dims = dims,
    cluster_k = cluster_k,
    cluster_methods = if (is.null(cluster_k)) character() else cluster_methods,
    hclust_method = hclust_method,
    spectral_sigma = spectral_sigma,
    spectral_nstart = spectral_nstart
  )
  class(out) <- c("MDistBenchmark", class(out))

  out
}

#' Extract pairwise distance benchmark comparisons
#'
#' Returns the pairwise diagnostics computed by [benchmark_mdist()].
#'
#' @param x An object returned by [benchmark_mdist()].
#'
#' @return A tibble with one row per unique pair of successful distance
#'   specifications. It contains pairwise MAD, symmetric relative distance,
#'   MDS congruence, alienation, and—when requested—one adjusted Rand index
#'   column per clustering method.
#'
#' @examples
#' \dontrun{
#' benchmark_comparisons(benchmark_result)
#' }
#'
#' @export
benchmark_comparisons <- function(x) {
  if (!inherits(x, "MDistBenchmark")) {
    stop("`x` must be an object returned by `benchmark_mdist()`.", call. = FALSE)
  }

  attr(x, "comparisons")
}

#' Plot pairwise distance benchmark comparisons
#'
#' Draws an annotated triangular heatmap of the pairwise diagnostics from
#' [benchmark_mdist()].
#'
#' @param object An object returned by [benchmark_mdist()].
#' @param metric Character string selecting `"mad"`, `"relative_distance"`,
#'   `"mds_congruence"`, `"alienation"`, or `"ari"`. A specific ARI column,
#'   such as `"ari_pam"`, can also be supplied.
#' @param cluster_method Optional clustering method used when `metric = "ari"`.
#'   If `NULL`, all available ARI metrics are shown, with one facet per
#'   clustering method.
#' @param digits Number of decimal places used for cell labels.
#' @param ... Currently unused.
#'
#' @return A `ggplot` object.
#'
#' @exportS3Method ggplot2::autoplot MDistBenchmark
autoplot.MDistBenchmark <- function(
    object,
    metric = c(
      "relative_distance", "mad", "alienation", "mds_congruence", "ari"
    ),
    cluster_method = NULL,
    digits = 2,
    ...
) {
  metric <- match.arg(
    metric,
    choices = c(
      "relative_distance", "mad", "alienation", "mds_congruence", "ari",
      "ari_pam", "ari_hclust", "ari_spectral"
    )
  )

  if (startsWith(metric, "ari_")) {
    cluster_method <- sub("^ari_", "", metric)
    metric <- "ari"
  }

  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      digits != as.integer(digits) || digits < 0L) {
    stop("`digits` must be a single non-negative integer.", call. = FALSE)
  }
  digits <- as.integer(digits)

  comparisons <- benchmark_comparisons(object)
  successful_labels <- attr(object, "method_labels")[object$ok]

  if (length(successful_labels) < 2L || nrow(comparisons) == 0L) {
    stop(
      "At least two successful distance specifications are required to plot ",
      "pairwise comparisons.",
      call. = FALSE
    )
  }

  ari_columns <- grep("^ari_", names(comparisons), value = TRUE)

  if (identical(metric, "ari")) {
    if (length(ari_columns) == 0L) {
      stop(
        "No ARI metrics are available. Re-run `benchmark_mdist()` with a ",
        "non-NULL `cluster_k`.",
        call. = FALSE
      )
    }

    if (!is.null(cluster_method)) {
      requested <- paste0("ari_", cluster_method)
      if (!requested %in% ari_columns) {
        stop(
          "ARI results are not available for clustering method `",
          cluster_method,
          "`.",
          call. = FALSE
        )
      }
      ari_columns <- requested
    }

    plot_data <- tidyr::pivot_longer(
      comparisons,
      cols = dplyr::all_of(ari_columns),
      names_to = "cluster_method",
      values_to = "value"
    ) |>
      dplyr::mutate(
        cluster_method = factor(
          sub("^ari_", "", .data$cluster_method),
          levels = sub("^ari_", "", ari_columns),
          labels = toupper(sub("^ari_", "", ari_columns))
        )
      )

    diagonal <- tidyr::crossing(
      method_1 = successful_labels,
      cluster_method = unique(plot_data$cluster_method)
    ) |>
      dplyr::transmute(
        method_1 = .data$method_1,
        method_2 = .data$method_1,
        cluster_method = .data$cluster_method,
        value = 1
      )

    plot_data <- dplyr::bind_rows(plot_data, diagonal)
    metric_label <- "Adjusted Rand index"
    facet_ari <- length(unique(plot_data$cluster_method)) > 1L
  } else {
    if (!metric %in% names(comparisons)) {
      stop("Metric `", metric, "` is not available.", call. = FALSE)
    }

    diagonal_value <- if (identical(metric, "mds_congruence")) 1 else 0
    plot_data <- comparisons |>
      dplyr::transmute(
        method_1 = .data$method_1,
        method_2 = .data$method_2,
        value = .data[[metric]]
      )
    diagonal <- tibble::tibble(
      method_1 = successful_labels,
      method_2 = successful_labels,
      value = diagonal_value
    )
    plot_data <- dplyr::bind_rows(plot_data, diagonal)

    metric_label <- switch(
      metric,
      mad = "Mean absolute difference",
      relative_distance = "Relative distance",
      mds_congruence = "MDS congruence",
      alienation = "Alienation coefficient"
    )
    facet_ari <- FALSE
  }

  plot_data <- plot_data |>
    dplyr::mutate(
      x_method = factor(.data$method_2, levels = successful_labels),
      y_method = factor(.data$method_1, levels = rev(successful_labels)),
      cell_label = ifelse(
        is.na(.data$value),
        "",
        sprintf(paste0("%.", digits, "f"), .data$value)
      )
    )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$x_method,
      y = .data$y_method,
      fill = .data$value
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$cell_label),
      size = 3
    ) +
    ggplot2::scale_fill_viridis_c(name = metric_label) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  if (facet_ari) {
    plot <- plot +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$cluster_method)
      )
  }

  plot
}
