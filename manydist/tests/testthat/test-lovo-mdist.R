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

lovo_plot_example <- function() {
  df <- tidyr::expand_grid(method = c("A", "B"),
                           variable = c("c_low", "c_high", "n_low", "n_high"))
  df$variable_type <- rep(c("categorical", "categorical", "numeric", "numeric"), 2)
  df$relative_distance <- c(0.1, 0.5, 0.8, 0.9, 0.8, 0.6, 0.1, 1)
  df$ari_pam <- df$relative_distance
  manydist:::MDistLOVOCompare$new(df, methods = list(A = list(), B = list()))
}

test_that("LOVO comparison reorders within type and preserves background bands", {
  result <- lovo_plot_example()
  p <- ggplot2::autoplot(result, reorder = TRUE)
  expect_identical(levels(p$data$variable), c("c_high", "c_low", "n_high", "n_low"))
  rects <- ggplot2::ggplot_build(p)$data[1:2]
  expect_equal(rects[[1]]$xmin, 0.5)
  expect_equal(rects[[1]]$xmax, 2.5)
  expect_equal(rects[[2]]$xmin, 2.5)
  expect_equal(rects[[2]]$xmax, 4.5)
  expect_identical(as.character(p$data$variable_type[p$data$variable %in%
                      levels(p$data$variable)[1:2]]), rep("categorical", 4))
  expect_identical(levels(ggplot2::autoplot(result)$data$variable),
                   c("c_low", "c_high", "n_low", "n_high"))
})

test_that("LOVO comparison ordering respects metric direction and global top_n", {
  result <- lovo_plot_example()
  p <- ggplot2::autoplot(result, metric = "ari_pam", reorder = TRUE)
  expect_identical(levels(p$data$variable), c("c_low", "c_high", "n_low", "n_high"))
  selected <- ggplot2::autoplot(result, reorder = TRUE, top_n = 3)
  expect_identical(levels(selected$data$variable), c("c_high", "c_low", "n_high"))
  only_numeric <- ggplot2::autoplot(result, reorder = TRUE, top_n = 1)
  expect_identical(levels(only_numeric$data$variable), "n_high")
  result$results$variable_type <- NULL
  global <- ggplot2::autoplot(result, reorder = TRUE)
  expect_identical(levels(global$data$variable), c("n_high", "c_high", "c_low", "n_low"))
})

test_that("single-method LOVO uses the same type-preserving ordering", {
  result <- lovo_mdist(lovo_example(), preset = "euclidean")
  result$results <- lovo_plot_example()$results |>
    dplyr::filter(method == "A")
  p <- ggplot2::autoplot(result, reorder = TRUE)
  expect_identical(levels(p$data$variable), c("c_high", "c_low", "n_high", "n_low"))
  ascending <- ggplot2::autoplot(result, metric = "ari_pam", reorder = TRUE)
  expect_identical(levels(ascending$data$variable), c("c_low", "c_high", "n_low", "n_high"))
})

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
