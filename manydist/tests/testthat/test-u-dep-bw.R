u_dep_bw_example <- function() {
  tibble::tibble(
    x1 = c(-2, -1, 0, 1, 2, 3),
    x2 = c(-3.0, -1.2, 0.4, 1.1, 3.2, 4.7),
    x3 = c(2.1, -0.5, 1.7, 3.4, -1.3, 0.8)
  )
}

manual_u_dep_bw <- function(x, ncomp = ncol(x), new_data = NULL) {
  x <- as.matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2, stats::sd)
  x_std <- sweep(sweep(x, 2, center, "-"), 2, scale, "/")
  pca <- stats::prcomp(x_std, center = FALSE, scale. = FALSE)

  retained <- seq_len(ncomp)
  train_scores <- sweep(
    pca$x[, retained, drop = FALSE],
    2,
    pca$sdev[retained],
    "/"
  )
  train_distance <- as.matrix(stats::dist(train_scores, method = "manhattan"))
  block_mean <- mean(train_distance[upper.tri(train_distance)])

  if (is.null(new_data)) {
    raw_distance <- train_distance
  } else {
    new_x <- as.matrix(new_data)
    new_x_std <- sweep(sweep(new_x, 2, center, "-"), 2, scale, "/")
    new_scores <- sweep(
      new_x_std %*% pca$rotation[, retained, drop = FALSE],
      2,
      pca$sdev[retained],
      "/"
    )

    raw_distance <- Reduce(
      `+`,
      Map(
        function(new_score, train_score) {
          abs(outer(new_score, train_score, "-"))
        },
        tibble::as_tibble(new_scores),
        tibble::as_tibble(train_scores)
      )
    )
  }

  raw_distance * ncol(x) / block_mean
}


test_that("u_dep_bw matches the block-wise whitened-PC definition", {
  x <- u_dep_bw_example()

  result <- mdist(x, preset = "u_dep_bw")
  observed <- as.matrix(result$distance)
  expected <- manual_u_dep_bw(x)

  expect_equal(observed, expected, tolerance = 1e-12)
  expect_equal(mean(observed[upper.tri(observed)]), ncol(x), tolerance = 1e-12)
  expect_equal(
    result$params$numeric_block_scale,
    ncol(x) / result$params$numeric_block_mean,
    tolerance = 1e-12
  )
  expect_identical(result$params$retained_ncomp, ncol(x))
})


test_that("u_dep_bw targets the original block size after dimension reduction", {
  x <- u_dep_bw_example()

  result <- mdist(x, preset = "u_dep_bw", ncomp = 2)
  observed <- as.matrix(result$distance)

  expect_equal(observed, manual_u_dep_bw(x, ncomp = 2), tolerance = 1e-12)
  expect_equal(mean(observed[upper.tri(observed)]), ncol(x), tolerance = 1e-12)
  expect_identical(result$params$retained_ncomp, 2L)
})


test_that("u_dep_bw holds the training block scale fixed for new data", {
  x <- u_dep_bw_example()
  new_data <- tibble::tibble(
    x1 = c(-0.5, 2.5),
    x2 = c(-0.2, 4.0),
    x3 = c(0.3, -0.7)
  )

  together <- mdist(x, new_data = new_data, preset = "u_dep_bw")
  separately <- mdist(x, new_data = new_data[1, ], preset = "u_dep_bw")

  expect_equal(
    as.matrix(together$distance),
    manual_u_dep_bw(x, new_data = new_data),
    tolerance = 1e-12
  )
  expect_equal(
    as.matrix(together$distance)[1, ],
    as.matrix(separately$distance)[1, ],
    tolerance = 1e-12
  )
})


test_that("adding u_dep_bw leaves the u_dep construction unchanged", {
  x <- u_dep_bw_example()

  preset_result <- mdist(x, preset = "u_dep")
  explicit_result <- mdist(
    x,
    preset = "custom",
    method_cat = "tvd",
    method_num = "pc_scores",
    commensurable = TRUE
  )

  expect_equal(
    as.matrix(preset_result$distance),
    as.matrix(explicit_result$distance),
    tolerance = 0
  )
})


test_that("u_dep_bw retains the u_dep categorical construction", {
  numerical <- u_dep_bw_example()
  categorical <- tibble::tibble(
    group_1 = factor(c("a", "a", "a", "b", "b", "b")),
    group_2 = factor(c("x", "x", "x", "y", "y", "y"))
  )

  mixed_result <- mdist(
    dplyr::bind_cols(numerical, categorical),
    preset = "u_dep_bw"
  )
  numerical_result <- mdist(numerical, preset = "u_dep_bw")
  categorical_result <- mdist(categorical, preset = "u_dep")

  expect_equal(
    as.matrix(mixed_result$distance),
    as.matrix(numerical_result$distance) +
      as.matrix(categorical_result$distance),
    tolerance = 1e-12
  )
})


test_that("u_dep_bw is registered as a response-aware preset", {
  registered <- dist_methods_tbl() |>
    dplyr::filter(.data$argument == "preset", .data$method == "u_dep_bw")

  expect_equal(nrow(registered), 1L)
  expect_true(registered$response_aware)
  expect_true(.step_mdist_is_response_aware("u_dep_bw", "matching"))
})


test_that("u_dep_bw rejects an undefined numerical block scale clearly", {
  constant_x <- tibble::tibble(x1 = 1:4, x2 = 2)

  expect_error(
    mdist(constant_x, preset = "u_dep_bw"),
    "positive finite standard deviations"
  )
})
