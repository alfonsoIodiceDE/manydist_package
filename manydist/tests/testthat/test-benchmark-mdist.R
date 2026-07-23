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
})

test_that("pairwise metrics are zero or one for identical specifications", {
  one_spec <- benchmark_specs("gower")
  specs <- dplyr::bind_rows(one_spec, one_spec)

  result <- benchmark_mdist(benchmark_example(), specs = specs)
  comparisons <- benchmark_comparisons(result)

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
  comparisons <- benchmark_comparisons(result)

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
    benchmark_comparisons(without_clusters)
  ))

  with_clusters <- benchmark_mdist(
    benchmark_example(),
    specs = specs,
    cluster_k = 3,
    cluster_methods = "pam"
  )
  comparisons <- benchmark_comparisons(with_clusters)

  expect_true("ari_pam" %in% names(comparisons))
  expect_equal(comparisons$ari_pam, 1)
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
})
