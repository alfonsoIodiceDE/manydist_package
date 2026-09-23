testthat::test_that("higher-dimensional generator is exactly paired", {
  reference <- moon_scaling_generate(
    seed = 26092100L, p_numeric = 2L, q_categorical = 2L,
    interaction_state = "signal"
  )
  signal <- moon_scaling_generate(
    seed = 26092100L, p_numeric = 10L, q_categorical = 5L,
    interaction_state = "signal"
  )
  null <- moon_scaling_generate(
    seed = 26092100L, p_numeric = 10L, q_categorical = 5L,
    interaction_state = "null"
  )

  testthat::expect_equal(dim(reference$data), c(360L, 4L))
  testthat::expect_equal(dim(signal$data), c(360L, 15L))
  testthat::expect_true(all(vapply(signal$data[1:10], is.numeric, logical(1))))
  testthat::expect_true(all(vapply(signal$data[11:15], is.factor, logical(1))))
  testthat::expect_false(anyNA(signal$data))
  testthat::expect_equal(
    as.vector(table(signal$latent$moon, signal$latent$band)),
    c(18L, 18L, 81L, 81L, 81L, 81L)
  )

  testthat::expect_identical(reference$data[1:2], signal$data[1:2])
  testthat::expect_identical(reference$data[3:4], signal$data[11:12])
  testthat::expect_identical(signal$data[1:10], null$data[1:10])
  testthat::expect_identical(signal$truth, null$truth)
  testthat::expect_identical(table(signal$data[11:15]), table(null$data[11:15]))
  testthat::expect_equal(signal$diagnostics$observed_band_agreement, 1)
  testthat::expect_lt(null$diagnostics$observed_band_agreement, 0.60)

  signal_delta <- getFromNamespace("cat_delta", "manydist")(
    signal$data[11:15], method_cat = "tvd"
  )$tvd
  null_delta <- getFromNamespace("cat_delta", "manydist")(
    null$data[11:15], method_cat = "tvd"
  )$tvd
  testthat::expect_equal(max(abs(signal_delta)), 0, tolerance = 1e-12)
  testthat::expect_equal(max(abs(null_delta)), 0, tolerance = 1e-12)
})

testthat::test_that("gated extension is nested in the pure package baseline", {
  signal <- moon_scaling_generate(
    seed = 26092100L, p_numeric = 2L, q_categorical = 2L,
    interaction_state = "signal"
  )
  null <- moon_scaling_generate(
    seed = 26092100L, p_numeric = 2L, q_categorical = 2L,
    interaction_state = "null"
  )
  fit_signal <- moon_scaling_gated_distance(
    signal$data, gamma = 1, prop_nn = 0.10,
    gate_permutations = 99L, gate_alpha = 0.01, gate_seed = 26092131L
  )
  fit_null <- moon_scaling_gated_distance(
    null$data, gamma = 1, prop_nn = 0.10,
    gate_permutations = 99L, gate_alpha = 0.01, gate_seed = 26092131L
  )
  direct <- manydist::mdist(
    signal$data,
    preset = "custom",
    method_num = "pc_scores",
    method_cat = "tvd",
    commensurable = FALSE,
    interaction = FALSE,
    ncomp = 2L
  )

  testthat::expect_equal(
    fit_signal$distance_no_interaction,
    as.matrix(direct$distance),
    tolerance = 0
  )
  testthat::expect_lt(fit_signal$diagnostics$decomposition_error, 1e-8)
  testthat::expect_true(fit_signal$gate$detected[["band"]])
  testthat::expect_false(fit_signal$gate$detected[["nuisance_1"]])
  testthat::expect_equal(fit_null$diagnostics$active_count, 0L)
  testthat::expect_identical(fit_null$distance, fit_null$distance_no_interaction)
  testthat::expect_lte(
    max(fit_signal$components$interaction_added),
    fit_signal$diagnostics$numeric_mean_distinct + 1e-10
  )
})
