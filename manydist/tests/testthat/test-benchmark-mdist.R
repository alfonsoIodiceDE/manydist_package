benchmark_example <- function() {
  tibble::tibble(
    group = factor(rep(c("a", "b", "c"), each = 4)),
    category = factor(rep(c("x", "y"), 6)),
    value_1 = c(1, 2, 1, 2, 5, 6, 5, 6, 9, 10, 9, 10),
    value_2 = c(2, 1, 2, 1, 6, 5, 6, 5, 10, 9, 10, 9)
  )
}

benchmark_specs <- function(presets = c("gower", "u_indep")) {
  all_dist_method_specs(
    mode = "presets_only",
    preset = presets
  )
}

test_that("specification grids use preset alone to distinguish custom settings", {
  specs <- all_dist_method_specs()
  expect_identical(names(specs), c("preset", "method_cat", "method_num", "commensurable"))
  custom <- dplyr::filter(specs, preset == "custom")
  expect_gt(nrow(custom), 0L)
  expect_false(anyNA(custom))
  presets <- all_dist_method_specs(mode = "presets_only")
  expect_false(any(presets$preset == "custom"))
  expect_true(all(is.na(presets$method_cat)))
})

test_that("benchmark dispatch is inferred from preset and ignores legacy type labels", {
  custom <- tibble::tibble(preset = "custom", method_cat = "matching",
                           method_num = "std", commensurable = TRUE)
  specs <- dplyr::bind_rows(benchmark_specs("gower"), custom)
  # Contradictory legacy labels cannot override the actual preset.
  specs$spec_type <- c("component", "preset")
  result <- benchmark_mdist(benchmark_example(), specs = specs)
  expect_true(all(result$ok))
  expect_false("spec_type" %in% names(result))
  expect_identical(result$result[[1]]$preset, "gower")
  expect_identical(result$result[[2]]$preset, "custom")
  expect_equal(result$result[[2]]$distance,
               mdist(benchmark_example(), method_cat = "matching",
                     method_num = "std", commensurable = TRUE)$distance)
  expect_match(attr(result, "method_labels")[2], "matching + std", fixed = TRUE)
  expect_output(print(result), "Methods:", fixed = TRUE)
  expect_error(benchmark_mdist(benchmark_example(),
                               specs = dplyr::mutate(custom, preset = NA_character_)),
               "nonmissing, nonempty")
})

test_that("benchmark_mdist retains its tibble interface", {
  result <- benchmark_mdist(
    benchmark_example(),
    specs = benchmark_specs()
  )

  expect_s3_class(result, "MDistBenchmark")
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("result", "ok", "error") %in% names(result)))
  expect_true(all(result$ok))

  selected <- dplyr::select(result, preset, ok)
  expect_equal(names(selected), c("preset", "ok"))
  expect_output(print(selected), "preset")

  reordered <- dplyr::arrange(result, dplyr::desc(.data$preset))
  expect_output(print(reordered), "preset")
})

test_that("benchmark printing is compact and reports failures", {
  successful <- benchmark_specs("gower") |>
    dplyr::mutate(label = "Gower")
  failed <- successful |>
    dplyr::mutate(
      label = "Invalid preset",
      preset = "not_a_preset"
    )

  result <- benchmark_mdist(
    benchmark_example(),
    specs = dplyr::bind_rows(successful, failed)
  )

  output <- capture.output(print(result))

  expect_true(any(grepl("MDistBenchmark", output, fixed = TRUE)))
  expect_true(any(grepl("specifications : 2", output, fixed = TRUE)))
  expect_true(any(grepl("successful     : 1", output, fixed = TRUE)))
  expect_true(any(grepl("failed         : 1", output, fixed = TRUE)))
  expect_true(any(grepl("Methods:", output, fixed = TRUE)))
  expect_true(any(grepl("Failures:", output, fixed = TRUE)))
  expect_true(any(grepl("Invalid preset", output, fixed = TRUE)))
  expect_true(any(grepl("summary(x)", output, fixed = TRUE)))
})

test_that("benchmark summary displays and returns the complete pairwise tibble", {
  result <- benchmark_mdist(
    benchmark_example(),
    specs = benchmark_specs(c("gower", "u_indep", "u_dep"))
  )
  stored <- attr(result, "comparisons")
  expect_output(
    visible_result <- withVisible(summary(result, n = 1L)),
    "Pairwise diagnostics",
    fixed = TRUE
  )
  expect_false(visible_result$visible)
  expect_s3_class(visible_result$value, "tbl_df")
  expect_identical(visible_result$value, stored)
  expect_equal(nrow(visible_result$value), 3L)
  expect_true(all(c("method_1", "method_2", "mad", "relative_distance",
                    "mds_congruence", "alienation") %in% names(stored)))
})

test_that("the obsolete benchmark extractor is not exported", {
  expect_false("benchmark_comparisons" %in% getNamespaceExports("manydist"))
})

test_that("benchmark summary handles fewer than two successful methods", {
  result <- benchmark_mdist(
    benchmark_example(),
    specs = benchmark_specs("gower")
  )

  metric_summary <- NULL
  expect_output(
    metric_summary <- summary(result),
    "No pairwise diagnostics are available",
    fixed = TRUE
  )
  expect_equal(nrow(metric_summary), 0L)
  expect_identical(metric_summary, attr(result, "comparisons"))
})

test_that("pairwise metrics are zero or one for identical specifications", {
  one_spec <- benchmark_specs("gower")
  specs <- dplyr::bind_rows(one_spec, one_spec)

  result <- benchmark_mdist(benchmark_example(), specs = specs)
  capture.output(comparisons <- summary(result))

  expect_equal(nrow(comparisons), 1L)
  expect_equal(comparisons$mad, 0)
  expect_equal(comparisons$relative_distance, 0)
  expect_equal(comparisons$mds_congruence, 1, tolerance = 1e-10)
  expect_equal(comparisons$alienation, 0, tolerance = 1e-7)
})

test_that("benchmark creates one row per unique method pair", {
  result <- benchmark_mdist(
    benchmark_example(),
    specs = benchmark_specs(c("gower", "u_indep", "u_dep"))
  )
  capture.output(comparisons <- summary(result))

  expect_equal(nrow(comparisons), choose(sum(result$ok), 2))
  expect_true(all(
    c("mad", "relative_distance", "mds_congruence", "alienation") %in%
      names(comparisons)
  ))
  expect_false(any(grepl("^ari_", names(comparisons))))
})

test_that("cluster_k controls optional pairwise ARI diagnostics", {
  one_spec <- benchmark_specs("gower")
  specs <- dplyr::bind_rows(one_spec, one_spec)

  without_clusters <- benchmark_mdist(
    benchmark_example(),
    specs = specs,
    cluster_k = NULL
  )
  expect_false("ari_pam" %in% names(
    attr(without_clusters, "comparisons")
  ))

  with_clusters <- benchmark_mdist(
    benchmark_example(),
    specs = specs,
    cluster_k = 3,
    cluster_methods = "pam"
  )
  capture.output(comparisons <- summary(with_clusters))

  expect_true("ari_pam" %in% names(comparisons))
  expect_equal(comparisons$ari_pam, 1)

  cluster_summary <- NULL
  expect_output(
    cluster_summary <- summary(with_clusters),
    "Pairwise diagnostics",
    fixed = TRUE
  )
  expect_true("ari_pam" %in% names(cluster_summary))
})

test_that("benchmark clustering arguments are validated", {
  expect_error(
    benchmark_mdist(
      benchmark_example(),
      specs = benchmark_specs(),
      cluster_k = 1
    ),
    "between 2 and nrow"
  )

  expect_error(
    benchmark_mdist(
      benchmark_example(),
      specs = benchmark_specs(),
      cluster_k = 3,
      cluster_methods = "unknown"
    ),
    "must be a subset"
  )
})

test_that("autoplot draws pairwise and clustering heatmaps", {
  one_spec <- benchmark_specs("gower")
  specs <- dplyr::bind_rows(one_spec, one_spec)

  result <- benchmark_mdist(
    benchmark_example(),
    specs = specs,
    cluster_k = 3,
    cluster_methods = "pam"
  )

  expect_s3_class(
    ggplot2::autoplot(result, metric = "relative_distance"),
    "ggplot"
  )
  expect_s3_class(
    ggplot2::autoplot(result, metric = "ari", cluster_method = "pam"),
    "ggplot"
  )

  for (metric in c(
    "relative_distance", "mad", "alienation", "mds_congruence", "ari", "ari_pam"
  )) {
    plot <- ggplot2::autoplot(result, metric = metric)
    fill_scale <- plot$scales$get_scales("fill")
    expect_equal(
      fill_scale$palette(c(0, 1)),
      c("#E76F51", "#008CFF"),
      info = metric
    )
    expect_equal(plot$layers[[2]]$aes_params$colour, "white", info = metric)
    expect_equal(plot$layers[[2]]$aes_params$fontface, "plain", info = metric)
    expect_false(any(plot$data$method_1 == plot$data$method_2), info = metric)
    expect_false(plot$scales$get_scales("x")$drop, info = metric)
    expect_false(plot$scales$get_scales("y")$drop, info = metric)
    if (startsWith(metric, "ari")) {
      expect_equal(fill_scale$limits, c(0, 1), info = metric)
    } else {
      expect_null(fill_scale$limits, info = metric)
    }
  }
})

test_that("ARI colours are fixed across observed ranges and retain negatives", {
  result <- benchmark_mdist(
    benchmark_example(),
    specs = benchmark_specs(),
    cluster_k = 3,
    cluster_methods = "pam"
  )

  capture.output(comparisons <- summary(result))
  comparisons$ari_pam <- 0.6
  attr(result, "comparisons") <- comparisons
  high_plot <- ggplot2::autoplot(result, metric = "ari")
  high_scale <- ggplot2::ggplot_build(high_plot)$plot$scales$get_scales("fill")
  expect_equal(high_scale$get_limits(), c(0, 1))
  expect_equal(high_scale$map(0.6), high_scale$palette(0.6))

  comparisons$ari_pam <- -0.2
  attr(result, "comparisons") <- comparisons
  negative_plot <- ggplot2::autoplot(result, metric = "ari")
  negative_build <- ggplot2::ggplot_build(negative_plot)
  negative_scale <- negative_build$plot$scales$get_scales("fill")
  expect_equal(negative_scale$map(0.6), high_scale$map(0.6))
  expect_equal(negative_build$data[[1]]$fill, "#E76F51")
  expect_equal(negative_build$data[[2]]$label, "-0.20")
  expect_equal(length(negative_build$layout$panel_params[[1]]$x$breaks), 2L)
  expect_equal(length(negative_build$layout$panel_params[[1]]$y$breaks), 2L)
})
