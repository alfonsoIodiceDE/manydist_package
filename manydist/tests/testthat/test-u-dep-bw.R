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

test_that("custom block_comm reproduces the preset and manual numerical definition", {
  x <- u_dep_bw_example()
  for (ncomp in c(2L, 3L)) {
    custom <- mdist(x, method_num = "pc_scores", commensurable = TRUE,
                    block_comm = TRUE, ncomp = ncomp)
    preset <- mdist(x, preset = "u_dep_bw", ncomp = ncomp)
    expect_equal(custom$distance, preset$distance, tolerance = 0)
    expect_equal(as.matrix(custom$distance), manual_u_dep_bw(x, ncomp = ncomp),
                 tolerance = 1e-12)
    expect_true(custom$params$block_comm)
    expect_true(preset$params$block_comm)
    expect_identical(preset$params$method_num, "pc_scores")
  }
  custom <- mdist(x, method_num = "pc_scores", block_comm = TRUE, threshold = 0.8)
  preset <- mdist(x, preset = "u_dep_bw", threshold = 0.8)
  expect_equal(custom$distance, preset$distance, tolerance = 0)
})

test_that("block_comm defaults to FALSE without changing componentwise PCs", {
  x <- u_dep_bw_example()
  default <- mdist(x, method_num = "pc_scores")
  explicit <- mdist(x, method_num = "pc_scores", block_comm = FALSE)
  expect_identical(default$distance, explicit$distance)
  expect_false(default$params$block_comm)
  expect_null(default$params$numeric_block_scale)
  expect_false(mdist(x, preset = "u_dep")$params$block_comm)
})

test_that("custom block scaling preserves categorical-method choice", {
  numeric <- u_dep_bw_example()
  categorical <- tibble::tibble(
    group_1 = factor(c("a", "a", "b", "b", "c", "c")),
    group_2 = factor(c("x", "y", "y", "x", "x", "y"))
  )
  mixed <- dplyr::bind_cols(numeric, categorical)
  for (method in c("matching", "tvd")) {
    full <- mdist(mixed, method_num = "pc_scores", method_cat = method,
                  block_comm = TRUE)
    numerical <- mdist(numeric, method_num = "pc_scores", block_comm = TRUE)
    cat <- mdist(categorical, method_cat = method)
    expect_equal(as.matrix(full$distance),
                 as.matrix(numerical$distance) + as.matrix(cat$distance),
                 tolerance = 1e-12)
  }
  expect_equal(
    mdist(mixed, method_num = "pc_scores", method_cat = "tvd",
          block_comm = TRUE)$distance,
    mdist(mixed, preset = "u_dep_bw")$distance,
    tolerance = 0
  )
})

test_that("custom block scaling fits numerical operators on training data only", {
  x <- u_dep_bw_example()
  test <- x[1:2, ] + 0.25
  train <- mdist(x, method_num = "pc_scores", block_comm = TRUE)
  together <- mdist(x, new_data = test, method_num = "pc_scores", block_comm = TRUE)
  alone <- mdist(x, new_data = test[1, ], method_num = "pc_scores", block_comm = TRUE)
  enlarged <- mdist(x, new_data = dplyr::bind_rows(test, test * 1000),
                    method_num = "pc_scores", block_comm = TRUE)
  expect_equal(as.matrix(together$distance), manual_u_dep_bw(x, new_data = test),
               tolerance = 1e-12)
  expect_equal(as.matrix(alone$distance)[1, ], as.matrix(together$distance)[1, ])
  expect_equal(as.matrix(enlarged$distance)[1:2, ], as.matrix(together$distance))
  expect_identical(train$params$numeric_block_scale, together$params$numeric_block_scale)
  reapplied <- mdist(x, new_data = x, method_num = "pc_scores", block_comm = TRUE)
  expect_equal(as.matrix(reapplied$distance), as.matrix(train$distance),
               tolerance = 1e-12)
})

test_that("block_comm is validated and non-custom presets retain control", {
  x <- u_dep_bw_example()
  for (value in list(NA, NULL, "TRUE", 1, c(TRUE, FALSE))) {
    expect_error(mdist(x, method_num = "pc_scores", block_comm = value),
                 "single TRUE or FALSE")
  }
  expect_error(mdist(x, method_num = "std", block_comm = TRUE),
               "requires.*pc_scores")
  expect_error(mdist(x, method_num = "pc_scores", commensurable = FALSE,
                     block_comm = TRUE), "commensurable = TRUE")
  expect_warning(gower <- mdist(x, preset = "gower", block_comm = TRUE),
                 "ignored.*block_comm")
  expect_false(gower$params$block_comm)
  single <- mdist(x[, 1, drop = FALSE], method_num = "pc_scores", block_comm = TRUE)
  expect_equal(mean(as.numeric(single$to_dist())), 1, tolerance = 1e-12)
})

test_that("step_mdist supports custom block scaling for pairwise and new data", {
  x <- u_dep_bw_example()
  test <- x[1:2, ] + 0.25
  rec <- recipes::recipe(~ ., data = x) |>
    step_mdist(recipes::all_predictors(), method_cat = "matching",
               method_num = "pc_scores", commensurable = TRUE, block_comm = TRUE) |>
    recipes::prep(training = x)
  expect_true(rec$steps[[1]]$block_comm)
  expect_true(rec$steps[[1]]$preprocessor$block_comm)
  baked <- recipes::bake(rec, new_data = test)
  expect_equal(unname(as.matrix(baked)),
               unname(manual_u_dep_bw(x, new_data = test)), tolerance = 1e-12)
  expect_equal(unname(as.matrix(recipes::bake(rec, new_data = test[1, ]))),
               unname(as.matrix(baked[1, ])), tolerance = 1e-12)
  pairwise <- recipes::recipe(~ ., data = x) |>
    step_mdist(recipes::all_predictors(), method_cat = "matching",
               method_num = "pc_scores", commensurable = TRUE,
               block_comm = TRUE, output = "pairwise") |>
    recipes::prep(training = x)
  expect_true(pairwise$steps[[1]]$block_comm)
  expect_equal(as.matrix(manydist::mdist(x, preset = "u_dep_bw")$distance),
               as.matrix(pairwise$steps[[1]]$preprocessor |>
                           .apply_mdist()), tolerance = 1e-12)
})

test_that("custom benchmark rows honor optional block_comm", {
  x <- u_dep_bw_example()
  specs <- tibble::tibble(preset = "custom", method_cat = "matching",
                          method_num = "pc_scores", commensurable = TRUE,
                          block_comm = c(FALSE, TRUE))
  result <- benchmark_mdist(x, specs = specs)
  expect_true(all(result$ok))
  expect_false(result$result[[1]]$params$block_comm)
  expect_true(result$result[[2]]$params$block_comm)
  expect_equal(result$result[[2]]$distance,
               mdist(x, method_cat = "matching", method_num = "pc_scores",
                     block_comm = TRUE)$distance, tolerance = 0)
  expect_match(attr(result, "method_labels")[2], "numerical block")
})
