knn_example <- function() {
  train <- data.frame(value = c(0, 2, 8, 10))
  test <- data.frame(value = c(1, 9))
  y <- factor(c("low", "low", "high", "high"), levels = c("low", "high"))
  list(train = train, test = test, y = y,
       d = mdist(train, new_data = test, preset = "gower"))
}

test_that("direct prediction accepts distances and tabular data", {
  ex <- knn_example()
  expected <- factor(c("low", "high"), levels = levels(ex$y))
  expect_identical(knn_dist(ex$d, ex$y, k = 1), expected)
  expect_identical(knn_dist(ex$d$distance, ex$y, k = 1), expected)
  expect_identical(knn_dist(as.data.frame(ex$d$distance), ex$y, k = 1), expected)
  expect_identical(knn_dist(ex$train, ex$y, new_data = ex$test, k = 1,
                           dist_fun = mdist, dist_args = list(preset = "gower")),
                   expected)
  engine <- fit_knn_dist(ex$train, ex$y, k = 1)
  expect_identical(knn_dist(ex$d, ex$y, k = 1),
                   predict_knn_dist_class(engine, ex$d$distance))
})

test_that("probabilities and regression use all neighbours", {
  ex <- knn_example()
  probabilities <- knn_dist(ex$d, ex$y, k = 2, type = "prob")
  expect_equal(probabilities, data.frame(low = c(1, 0), high = c(0, 1)))
  expect_equal(rowSums(probabilities), c(1, 1))
  expect_equal(knn_dist(ex$d, ex$train$value, k = 2, type = "numeric"), c(1, 9))
  expect_equal(knn_dist(ex$train, ex$train$value, new_data = ex$test, k = 2,
                        dist_fun = mdist, dist_args = list(preset = "gower"),
                        type = "numeric"), c(1, 9))
  expect_equal(knn_dist(ex$train, ex$y, new_data = ex$test, k = 2,
                        dist_fun = mdist, dist_args = list(preset = "gower"),
                        type = "prob"), probabilities)
})

test_that("single test rows and single outcome levels retain the right shape", {
  D <- matrix(c(0, 1, 2), nrow = 1)
  y <- factor(c("a", "b", "b"))
  expect_equal(as.character(knn_dist(D, y, k = 3)), "b")
  expect_equal(knn_dist(D, y, k = 3, type = "prob"),
               data.frame(a = 1/3, b = 2/3))
  expect_equal(knn_dist(D, c(1, 4, 10), k = 3, type = "numeric"), 5)
  one_class <- factor(rep("a", 3))
  expect_equal(knn_dist(rbind(D, D), one_class, k = 2, type = "prob"),
               data.frame(a = c(1, 1)))
})

test_that("empty test matrices produce empty predictions", {
  D <- matrix(numeric(), nrow = 0, ncol = 2)
  y <- factor(c("a", "b"))
  expect_identical(knn_dist(D, y, k = 1), factor(character(), levels = levels(y)))
  expect_equal(dim(knn_dist(D, y, k = 1, type = "prob")), c(0L, 2L))
  expect_identical(knn_dist(D, c(1, 2), k = 1, type = "numeric"), numeric())
})

test_that("ties follow training order and factor level order", {
  y <- factor(c("b", "a"), levels = c("a", "b"))
  D <- matrix(c(1, 1), nrow = 1)
  expect_equal(as.character(knn_dist(D, y, k = 1)), "b")
  expect_equal(as.character(knn_dist(D, y, k = 2)), "a")
})

test_that("invalid direct inputs give informative errors", {
  ex <- knn_example()
  for (k in list(0, -1, 1.5, 5, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(knn_dist(ex$d, ex$y, k = k), "integer between 1")
  }
  expect_error(knn_dist(ex$d, ex$y[-1], k = 1), "one column per")
  expect_error(knn_dist(ex$d, as.character(ex$y), k = 1), "factor outcomes")
  expect_error(knn_dist(ex$d, ex$y, k = 1, type = "numeric"), "numeric outcomes")
  missing_y <- ex$y
  missing_y[1] <- NA
  expect_error(knn_dist(ex$d, missing_y, k = 1), "nonmissing")
  for (bad in c(NA_real_, Inf, -1)) {
    D <- as.matrix(ex$d$distance)
    D[1, 1] <- bad
    expect_error(knn_dist(D, ex$y, k = 1), "finite, and nonnegative")
  }
  expect_error(knn_dist(matrix("a", nrow = 2, ncol = 4), ex$y, k = 1),
               "numeric, finite")
  expect_error(knn_dist(1:4, ex$y, k = 1), "test-to-training matrix")
  expect_error(knn_dist(ex$d, ex$y, new_data = ex$test, k = 1),
               "leave `new_data`")
  expect_error(knn_dist(ex$train, ex$y, k = 1, dist_fun = "mdist"),
               "function or NULL")
  expect_error(knn_dist(ex$train, ex$y, k = 1, dist_fun = mdist),
               "tabular test predictors")
  expect_error(knn_dist(ex$train[-1, , drop = FALSE], ex$y, new_data = ex$test,
                        k = 1, dist_fun = mdist), "one row per outcome")
  expect_error(knn_dist(ex$train, ex$y, new_data = ex$test, k = 1,
                        dist_fun = mdist, dist_args = list("gower")), "named list")
  expect_error(knn_dist(ex$train, ex$y, new_data = ex$test, k = 1,
                        dist_fun = mdist, dist_args = list(new_data = ex$test)),
               "not `dist_args`")
  wrong_rows <- function(x, new_data) matrix(0, nrow = 3, ncol = nrow(x))
  expect_error(knn_dist(ex$train, ex$y, new_data = ex$test, k = 1,
                        dist_fun = wrong_rows), "one row per test")
})

test_that("tabular prediction calls the distance function with separate datasets", {
  ex <- knn_example()
  called <- 0L
  distance_function <- function(x, new_data, scale) {
    called <<- called + 1L
    expect_identical(x, ex$train)
    expect_identical(new_data, ex$test)
    as.data.frame(scale * abs(outer(new_data$value, x$value, "-")))
  }
  expect_identical(knn_dist(ex$train, ex$y, new_data = ex$test, k = 1,
                           dist_fun = distance_function, dist_args = list(scale = 2)),
                   knn_dist(ex$d, ex$y, k = 1))
  expect_equal(called, 1L)
})

test_that("test extremes do not change training-only Gower scaling", {
  ex <- knn_example()
  extended <- data.frame(value = c(ex$test$value, 100))
  original <- knn_dist(ex$train, ex$y, new_data = ex$test, k = 2,
                        dist_fun = mdist, dist_args = list(preset = "gower"),
                        type = "prob")
  extra <- knn_dist(ex$train, ex$y, new_data = extended, k = 2,
                     dist_fun = mdist, dist_args = list(preset = "gower"),
                     type = "prob")
  expect_equal(extra[1:2, ], original)
})

test_that("direct and tidymodels workflow predictions agree", {
  ex <- knn_example()
  data <- ex$train
  data$outcome <- ex$y
  rec <- recipes::recipe(outcome ~ ., data = data) |>
    step_mdist(recipes::all_predictors(), preset = "gower",
               output = "distance_to_training")
  spec <- nearest_neighbor_dist(mode = "classification", neighbors = 2) |>
    parsnip::set_engine("manydist")
  wf <- workflows::workflow() |>
    workflows::add_recipe(rec) |>
    workflows::add_model(spec)
  fit <- workflows::fit(wf, data = data)
  expect_identical(predict(fit, new_data = ex$test)$.pred_class,
                   knn_dist(ex$d, ex$y, k = 2))
})

test_that("response selection supports names and excludes test outcomes", {
  ex <- knn_example()
  train <- ex$train
  train$outcome <- ex$y
  test <- ex$test
  test$outcome <- factor(c("high", "low"), levels = levels(ex$y))
  expected <- knn_dist(ex$d, ex$y, k = 2)
  from_name <- knn_dist(train, response = "outcome", new_data = test, k = 2,
                         dist_fun = mdist, dist_args = list(preset = "gower"))
  from_symbol <- knn_dist(train, response = outcome, new_data = ex$test, k = 2,
                           dist_fun = mdist, dist_args = list(preset = "gower"))
  expect_identical(from_name, expected)
  expect_identical(from_symbol, expected)
  test$outcome <- NA
  expect_identical(knn_dist(train, response = "outcome", new_data = test, k = 2,
                            dist_fun = mdist, dist_args = list(preset = "gower")),
                   expected)
  checked_distance <- function(x, new_data) {
    expect_identical(x, ex$train)
    expect_identical(new_data, ex$test)
    mdist(x, new_data = new_data, preset = "gower")
  }
  expect_identical(knn_dist(train, response = "outcome", new_data = test,
                            k = 2, dist_fun = checked_distance), expected)
  expect_error(knn_dist(train, ex$y, response = "outcome", new_data = test,
                        k = 2, dist_fun = mdist), "either `y` or `response`")
  expect_error(knn_dist(ex$d, response = "outcome", k = 2),
               "tabular training data")
  expect_error(knn_dist(train, response = c(value, outcome), new_data = test,
                        k = 2, dist_fun = mdist), "exactly one outcome")
  expect_error(knn_dist(ex$d, k = 2), "nonempty vector")
})

test_that("response selection works for numeric regression", {
  ex <- knn_example()
  train <- ex$train
  train$outcome <- c(0, 2, 8, 10)
  expect_equal(knn_dist(train, response = "outcome", new_data = ex$test, k = 2,
                        dist_fun = mdist, dist_args = list(preset = "gower"),
                        type = "numeric"), c(1, 9))
})
