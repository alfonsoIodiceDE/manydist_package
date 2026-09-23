lovo_example <- function() {
  tibble::tibble(
    x = c(1, 2, 3, 4, 5, 6, 7, 8),
    y = c(1, 4, 2, 8, 5, 7, 3, 6),
    z = c(8, 1, 6, 2, 7, 3, 5, 4)
  )
}

lovo_methods <- function() {
  list(
    euclidean_std = list(preset = "euclidean"),
    euclidean_range = list(preset = "euclidean", method_num = "range")
  )
}

test_that("LOVO skips MDS diagnostics by default", {
  result <- lovo_mdist(lovo_example(), preset = "euclidean")

  expect_false(result$mds)
  expect_null(result$dims)
  expect_null(result$base_mds)
  expect_false(any(c(
    "cc_importance", "mds_congruence", "ac_importance"
  ) %in% names(result$results)))
  expect_output(print(result), "MDS diagnostics : not computed", fixed = TRUE)
  expect_error(
    ggplot2::autoplot(result, metric = "mds_congruence"),
    "Re-run with `mds = TRUE`",
    fixed = TRUE
  )
})

test_that("mds requests MDS diagnostics with the selected dimensions", {
  result <- lovo_mdist(
    lovo_example(),
    preset = "euclidean",
    mds = TRUE,
    dims = 2
  )

  expect_true(result$mds)
  expect_equal(result$dims, 2L)
  expect_equal(dim(result$base_mds), c(nrow(lovo_example()), 2L))
  expect_true(all(c(
    "cc_importance", "mds_congruence", "ac_importance"
  ) %in% names(result$results)))
  expect_s3_class(
    ggplot2::autoplot(result, metric = "mds_congruence"),
    "ggplot"
  )
})

test_that("LOVO summary omits response and MDS invocation details", {
  result <- lovo_mdist(
    lovo_example(),
    preset = "gower",
    mds = TRUE,
    dims = 2
  )

  output <- capture.output(summary(result))

  expect_false(any(grepl("response used", output, fixed = TRUE)))
  expect_false(any(grepl("MDS diagnostics", output, fixed = TRUE)))
  expect_false(any(grepl("dims", output, fixed = TRUE)))
  expect_true(any(grepl("preset : gower", output, fixed = TRUE)))
})

test_that("LOVO validates the MDS controls", {
  expect_error(
    lovo_mdist(lovo_example(), preset = "euclidean", mds = NA),
    "mds must be a single TRUE/FALSE value",
    fixed = TRUE
  )
  expect_error(
    lovo_mdist(
      lovo_example(),
      preset = "euclidean",
      mds = TRUE,
      dims = 0
    ),
    "dims must be a single positive integer",
    fixed = TRUE
  )

  expect_warning(
    result <- lovo_mdist(
      lovo_example(),
      preset = "euclidean",
      dims = 1
    ),
    "`dims` is ignored when `mds = FALSE`",
    fixed = TRUE
  )
  expect_false(result$mds)
})

test_that("LOVO comparison supports default-off and requested MDS", {
  without_mds <- compare_lovo_mdist(
    lovo_example(),
    methods = lovo_methods()
  )

  expect_false(without_mds$mds)
  expect_null(without_mds$dims)
  expect_false("mds_congruence" %in% names(without_mds$results))

  without_summary <- NULL
  expect_output(
    without_summary <- summary(without_mds),
    "MDS diagnostics : not computed",
    fixed = TRUE
  )
  expect_false(any(c("mds_min", "mds_max", "mds_mean") %in%
                     names(without_summary)))

  with_mds <- compare_lovo_mdist(
    lovo_example(),
    methods = lovo_methods(),
    mds = TRUE,
    dims = 2
  )

  expect_true(with_mds$mds)
  expect_equal(with_mds$dims, 2L)
  expect_true("mds_congruence" %in% names(with_mds$results))

  with_summary <- NULL
  expect_output(
    with_summary <- summary(with_mds),
    "MDS diagnostics : 2 dimensions",
    fixed = TRUE
  )
  expect_true(all(c("mds_min", "mds_max", "mds_mean") %in%
                    names(with_summary)))
})

test_that("LOVO comparison warns when dims is supplied without MDS", {
  expect_warning(
    result <- compare_lovo_mdist(
      lovo_example(),
      methods = lovo_methods(),
      dims = 1
    ),
    "`dims` is ignored when `mds = FALSE`",
    fixed = TRUE
  )
  expect_false(result$mds)
})

test_that("MDS and clustering diagnostics can be requested independently", {
  skip_if_not_installed("mclust")

  clusters_only <- lovo_mdist(
    lovo_example(),
    preset = "euclidean",
    cluster_k = 2,
    cluster_methods = "pam"
  )

  expect_true("pam_importance" %in% names(clusters_only$results))
  expect_false("mds_congruence" %in% names(clusters_only$results))

  both <- lovo_mdist(
    lovo_example(),
    preset = "euclidean",
    mds = TRUE,
    dims = 2,
    cluster_k = 2,
    cluster_methods = "pam"
  )

  expect_true(all(c("pam_importance", "mds_congruence") %in%
                    names(both$results)))
})
