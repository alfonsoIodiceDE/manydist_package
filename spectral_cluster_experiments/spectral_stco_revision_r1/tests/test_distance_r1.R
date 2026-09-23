testthat::test_that("the numerical block is exactly the pinned manydist u_dep_bw block", {
  data <- data.frame(
    x1 = seq(-2, 2, length.out = 12),
    x2 = c(-1.9, -1.1, -0.8, -0.2, 0.4, 0.9, 1.4, 1.7, 2.2, 2.9, 3.1, 3.8),
    c1 = factor(rep(c("a", "b"), each = 6)),
    c2 = factor(rep(c("u", "v", "u"), each = 4))
  )

  result <- scmix_distance_r1(
    data,
    interaction_gate = "none",
    strict_provenance = TRUE
  )
  expected <- manydist::mdist(data[c("x1", "x2")], preset = "u_dep_bw")

  testthat::expect_equal(
    result$components$numeric$scaled,
    as.matrix(expected$distance),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  testthat::expect_equal(
    mean(result$components$numeric$scaled[upper.tri(result$components$numeric$scaled)]),
    2,
    tolerance = 1e-12
  )
})

testthat::test_that("the final categorical blocks use raw bounded IFCS deltas", {
  set.seed(20260918)
  n <- 48
  c1 <- factor(rep(c("a", "b"), each = n / 2))
  c2 <- factor(ifelse(c1 == "a", "u", "v"))
  data <- data.frame(
    x1 = as.numeric(c1) + stats::rnorm(n, sd = 0.15),
    x2 = as.numeric(c1) + stats::rnorm(n, sd = 0.20),
    c1 = c1,
    c2 = c2
  )

  result <- scmix_distance_r1(data, interaction_gate = "none")
  testthat::expect_equal(result$config$rho, 1 / 2)

  for (component in result$components$categorical) {
    testthat::expect_true(component$tvd_available)
    testthat::expect_true(component$interaction_available)
    testthat::expect_gte(min(component$tvd_raw), 0)
    testthat::expect_lte(max(component$tvd_raw), 1)
    testthat::expect_gte(min(component$interaction_raw), 0)
    testthat::expect_lte(max(component$interaction_raw), 1)
    testthat::expect_equal(
      component$combined,
      component$tvd_raw + 0.5 * component$interaction_raw,
      tolerance = 1e-12
    )
    testthat::expect_equal(
      component$effective_weights,
      c(tvd = 1, interaction = 0.5),
      tolerance = 0
    )

    # Retained only for diagnostics; these matrices do not enter the final D.
    testthat::expect_equal(
      mean(component$tvd_scaled[upper.tri(component$tvd_scaled)]),
      1,
      tolerance = 1e-12
    )
    testthat::expect_equal(
      mean(component$interaction_scaled[upper.tri(component$interaction_scaled)]),
      1,
      tolerance = 1e-12
    )
  }
})

testthat::test_that("a single categorical variable follows the same fixed formula", {
  theta <- seq(0, pi, length.out = 30)
  data <- rbind(
    data.frame(x1 = cos(theta), x2 = sin(theta), group = factor("upper", levels = c("upper", "lower"))),
    data.frame(x1 = 1 - cos(theta), x2 = 0.45 - sin(theta), group = factor("lower", levels = c("upper", "lower")))
  )

  result <- scmix_distance_r1(
    data,
    prop_nn = 0.10,
    interaction_gate = "none"
  )
  component <- result$components$categorical$group

  testthat::expect_equal(result$config$rho, 1)
  testthat::expect_false(component$tvd_available)
  testthat::expect_true(component$interaction_available)
  testthat::expect_identical(
    component$availability_reason,
    "interaction_detected_tvd_zero"
  )
  testthat::expect_equal(component$effective_weights, c(tvd = 1, interaction = 1))
  testthat::expect_equal(component$combined, component$interaction_raw)
  testthat::expect_equal(
    result$distance_no_interaction,
    result$components$numeric$scaled,
    tolerance = 1e-12
  )
})

testthat::test_that("a zero interaction keeps its fixed coefficient at zero", {
  n <- 20
  alternating <- factor(rep(c("a", "b"), n / 2))
  data <- data.frame(
    x1 = seq_len(n),
    c1 = alternating,
    c2 = factor(ifelse(alternating == "a", "u", "v"))
  )

  result <- scmix_distance_r1(
    data,
    prop_nn = 0.10,
    interaction_gate = "none"
  )
  component <- result$components$categorical$c1

  testthat::expect_true(component$tvd_available)
  testthat::expect_false(component$interaction_available)
  testthat::expect_identical(
    component$availability_reason,
    "interaction_empirically_zero"
  )
  testthat::expect_equal(component$effective_weights, c(tvd = 1, interaction = 0.5))
  testthat::expect_equal(component$combined, component$tvd_raw)
})

testthat::test_that("the final and no-interaction distances are exact raw compositions", {
  data <- data.frame(
    x1 = c(-3, -2, -1, 0, 1, 2, 3, 4),
    x2 = c(-2.5, -1.7, -1.2, 0.3, 0.6, 2.4, 2.7, 4.1),
    c1 = factor(c("a", "a", "a", "b", "b", "b", "a", "b")),
    c2 = factor(c("u", "u", "v", "v", "v", "u", "u", "v"))
  )

  result <- scmix_distance_r1(data, interaction_gate = "none")
  expected <- result$components$numeric$scaled + result$components$categorical_sum
  expected_no_interaction <- result$components$numeric$scaled +
    Reduce(
      `+`,
      lapply(result$components$categorical, `[[`, "tvd_raw"),
      init = matrix(0, nrow(data), nrow(data))
    )

  testthat::expect_equal(result$distance, expected, tolerance = 1e-12)
  testthat::expect_equal(
    result$distance_no_interaction,
    expected_no_interaction,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  testthat::expect_equal(result$distance, t(result$distance), tolerance = 1e-12)
  testthat::expect_equal(unname(diag(result$distance)), rep(0, nrow(data)), tolerance = 0)
  testthat::expect_identical(result$provenance$package_version, "0.5.2")
  testthat::expect_identical(
    result$provenance$source_commit,
    "ab27241f1a0846bb13087f804e23fd52746d1578"
  )
})

testthat::test_that("rho defaults to 1/q and permits pre-specified sensitivity values", {
  data <- data.frame(
    x = seq_len(12),
    c1 = factor(rep(c("a", "b"), each = 6)),
    c2 = factor(rep(c("u", "v"), 6))
  )

  default_fit <- scmix_distance_r1(data, interaction_gate = "none")
  sensitivity_fit <- scmix_distance_r1(
    data,
    rho = 0.25,
    interaction_gate = "none"
  )

  testthat::expect_equal(default_fit$config$rho, 0.5)
  testthat::expect_true(default_fit$config$rho_is_default)
  testthat::expect_equal(sensitivity_fit$config$rho, 0.25)
  testthat::expect_false(sensitivity_fit$config$rho_is_default)
  testthat::expect_error(
    scmix_distance_r1(data, rho = -0.01, interaction_gate = "none"),
    "in [0, 1]",
    fixed = TRUE
  )
  testthat::expect_error(
    scmix_distance_r1(data, rho = 1.01, interaction_gate = "none"),
    "in [0, 1]",
    fixed = TRUE
  )
})

testthat::test_that("the permutation gate nests the no-interaction distance", {
  n_per_moon <- 45L
  theta <- ((seq_len(n_per_moon) - 0.5) / n_per_moon) * pi
  continuous <- rbind(
    data.frame(x1 = cos(theta), x2 = sin(theta), moon = "upper"),
    data.frame(x1 = 1 - cos(theta), x2 = 0.5 - sin(theta), moon = "lower")
  )
  continuous$moon <- factor(continuous$moon, levels = c("upper", "lower"))
  thirds <- rep(c("A", "B", "C"), each = n_per_moon / 3L)

  signal <- continuous
  signal$band <- factor(rep(thirds, 2L), levels = c("A", "B", "C"))
  signal$nuisance <- factor(rep(c("N1", "N2", "N3"), length.out = nrow(signal)))

  null <- continuous
  set.seed(11)
  null$band <- factor(c(sample(thirds), sample(thirds)), levels = c("A", "B", "C"))
  set.seed(12)
  null$nuisance <- factor(sample(rep(c("N1", "N2", "N3"), length.out = nrow(null))))

  signal_fit <- scmix_distance_r1(
    signal[c("x1", "x2", "band", "nuisance")],
    gate_permutations = 19L,
    gate_alpha = 0.05,
    gate_seed = 31L
  )
  null_fit <- scmix_distance_r1(
    null[c("x1", "x2", "band", "nuisance")],
    gate_permutations = 19L,
    gate_alpha = 0.05,
    gate_seed = 32L
  )

  testthat::expect_true(signal_fit$diagnostics$interaction_detected[["band"]])
  testthat::expect_false(any(null_fit$diagnostics$interaction_detected))
  testthat::expect_equal(
    null_fit$distance,
    null_fit$distance_no_interaction,
    tolerance = 0
  )
  testthat::expect_lte(
    max(
      signal_fit$distance - signal_fit$distance_no_interaction
    ),
    1 + 1e-12
  )
})
