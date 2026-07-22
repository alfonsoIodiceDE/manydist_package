response_example <- function() {
  tibble::tibble(
    outcome = factor(c(
      "low", "low", "low", "low",
      "low", "low", "high", "high",
      "high", "high", "high", "high"
    )),
    region = factor(rep(c("a", "b", "c"), each = 4)),
    score = c(1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6)
  )
}

test_that("step_mdist uses a single recipe outcome for response-aware presets", {
  dat <- response_example()

  rec <- recipes::recipe(outcome ~ ., data = dat) |>
    step_mdist(
      recipes::all_predictors(),
      preset = "u_mix",
      output = "distance_to_training"
    )

  rec_prep <- recipes::prep(rec, training = dat)
  mdist_step <- rec_prep$steps[[1]]

  expect_identical(mdist_step$response_col, "outcome")
  expect_identical(mdist_step$preprocessor$response_col, "outcome")
  expect_identical(mdist_step$preprocessor$y, dat$outcome)
  expect_identical(mdist_step$preprocessor$method_cat, "tvd")
  expect_false("outcome" %in% names(mdist_step$preprocessor$x_train))
})

test_that("baking new data never requires or uses its outcome", {
  dat <- response_example()

  rec_prep <- recipes::recipe(outcome ~ ., data = dat) |>
    step_mdist(recipes::all_predictors(), preset = "u_mix") |>
    recipes::prep(training = dat)

  new_predictors <- dat[1:3, c("region", "score")]
  new_with_outcome <- dat[1:3, ]
  new_with_outcome$outcome <- factor(
    rev(new_with_outcome$outcome),
    levels = levels(dat$outcome)
  )

  baked_without <- recipes::bake(rec_prep, new_data = new_predictors)
  baked_with <- recipes::bake(rec_prep, new_data = new_with_outcome)

  distance_names <- grep("^dist_", names(baked_without), value = TRUE)
  expect_equal(
    baked_without[, distance_names],
    baked_with[, distance_names]
  )
})

test_that("response use can be disabled explicitly", {
  dat <- response_example()

  rec <- recipes::recipe(outcome ~ ., data = dat) |>
    step_mdist(
      recipes::all_predictors(),
      preset = "u_mix",
      response_used = FALSE
    )

  expect_warning(
    rec_prep <- recipes::prep(rec, training = dat),
    "Without a response, 'tvd' requires >1 categorical variable"
  )
  mdist_step <- rec_prep$steps[[1]]

  expect_null(mdist_step$response_col)
  expect_null(mdist_step$preprocessor$y)
  expect_identical(mdist_step$preprocessor$method_cat, "matching")
})

test_that("non-response-aware presets do not consume the recipe outcome", {
  dat <- response_example()

  rec_prep <- recipes::recipe(outcome ~ ., data = dat) |>
    step_mdist(recipes::all_predictors(), preset = "gower") |>
    recipes::prep(training = dat)

  mdist_step <- rec_prep$steps[[1]]

  expect_null(mdist_step$response_col)
  expect_null(mdist_step$preprocessor$y)
})
