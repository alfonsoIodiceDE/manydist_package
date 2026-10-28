test_that("non-PC numerical preprocessing omits PCA print fields", {
  x <- iris[, 1:4]
  for (method in c("none", "std", "range", "robust")) {
    object <- mdist(x, method_num = method, commensurable = FALSE,
                    ncomp = 2, threshold = 0.8)
    output <- capture.output(returned <- object$print())
    expect_identical(returned, object)
    expect_false(any(grepl("principal components|inertia threshold", output)))
    selected <- capture.output(object$print(show = c("method_num", "ncomp", "threshold")))
    expect_false(any(grepl("principal components|inertia threshold", selected)))
    expect_true(any(grepl("numerical preprocessing", output, fixed = TRUE)))
  }
})

test_that("PC scores show component counts but omit NULL thresholds", {
  object <- mdist(iris[, 1:4], method_num = "pc_scores", ncomp = 2)
  output <- capture.output(print(object))
  expect_true(any(grepl("number of principal components: 2", output, fixed = TRUE)))
  expect_false(any(grepl("inertia threshold", output, fixed = TRUE)))
  selected <- capture.output(object$print(show = c("ncomp", "threshold")))
  expect_false(any(grepl("inertia threshold", selected, fixed = TRUE)))
})

test_that("PC scores show non-NULL thresholds", {
  object <- mdist(iris[, 1:4], method_num = "pc_scores", threshold = 0.8)
  output <- capture.output(print(object))
  expect_true(any(grepl("inertia threshold: 0.8", output, fixed = TRUE)))
  expect_identical(object$params$threshold, 0.8)
})

test_that("blockwise PCs retain their relevant count and scale fields", {
  object <- mdist(iris[, 1:4], preset = "u_dep_bw", ncomp = 2)
  output <- capture.output(print(object))
  expect_true(any(grepl("retained principal components: 2", output, fixed = TRUE)))
  expect_true(any(grepl("numerical block scale factor", output, fixed = TRUE)))
  expect_false(any(grepl("inertia threshold", output, fixed = TRUE)))
})
