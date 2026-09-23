testthat::test_that("c = 1 reproduces the frozen Gaussian affinity", {
  distance <- as.matrix(stats::dist(matrix(c(0, 0, 1, 0, 0, 2), ncol = 2, byrow = TRUE)))
  expected <- full_affinity(distance, "gaussian", list())
  observed <- moon_bandwidth_gaussian_affinity(distance, 1)

  testthat::expect_identical(observed$matrix, expected$matrix)
  testthat::expect_identical(observed$scale, expected$scale)
  testthat::expect_equal(observed$base_scale, stats::median(distance[upper.tri(distance)]))
})

testthat::test_that("the multiplier changes only the declared global scale", {
  distance <- as.matrix(stats::dist(matrix(c(0, 0, 1, 0, 0, 2), ncol = 2, byrow = TRUE)))
  narrow <- moon_bandwidth_gaussian_affinity(distance, 0.5)
  broad <- moon_bandwidth_gaussian_affinity(distance, 2)

  testthat::expect_equal(narrow$scale, 0.5 * narrow$base_scale)
  testthat::expect_equal(broad$scale, 2 * broad$base_scale)
  testthat::expect_true(all(narrow$matrix[upper.tri(narrow$matrix)] <=
                            broad$matrix[upper.tri(broad$matrix)]))
})

testthat::test_that("invalid multipliers are rejected", {
  distance <- as.matrix(stats::dist(matrix(1:6, ncol = 2)))
  testthat::expect_error(moon_bandwidth_gaussian_affinity(distance, 0), "positive")
  testthat::expect_error(moon_bandwidth_gaussian_affinity(distance, NA_real_), "positive")
})

