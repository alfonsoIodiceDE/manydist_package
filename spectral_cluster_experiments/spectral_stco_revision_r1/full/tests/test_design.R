testthat::test_that("the frozen design expands to the declared task count", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  tasks <- full_expand_tasks(design)

  testthat::expect_length(tasks, 1650L)
  testthat::expect_equal(
    length(tasks),
    as.integer(design$execution$expected_task_count)
  )
})

testthat::test_that("the nested interaction gate is fully pre-specified", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  gate <- design$distance$categorical$interaction_gate

  testthat::expect_identical(
    design$distance$categorical$composition,
    "package_baseline_plus_gamma_times_numeric_mean_times_active_interaction_average"
  )
  testthat::expect_identical(
    design$distance$categorical$active_component_rule,
    "mean_over_multiplicity_gated_active_variables"
  )
  testthat::expect_equal(design$distance$defaults$gamma, 1)
  testthat::expect_identical(gate$method, "permutation_max")
  testthat::expect_equal(gate$permutations, 99)
  testthat::expect_equal(gate$alpha, 0.01)
  testthat::expect_equal(gate$seed_offset, 31)
})

testthat::test_that("the full runner applies the active-average construction", {
  numeric_distance <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3L, 3L)
  categorical_distance <- matrix(c(0, 0.2, 0.4, 0.2, 0, 0.6, 0.4, 0.6, 0), 3L, 3L)
  active_average <- matrix(c(0, 0.5, 1, 0.5, 0, 0.25, 1, 0.25, 0), 3L, 3L)
  level_delta <- matrix(c(0, 0.7, 0.7, 0), 2L, 2L)
  observation_delta <- matrix(c(0, 0.5, 1, 0.5, 0, 0.25, 1, 0.25, 0), 3L, 3L)
  fit <- list(
    components = list(
      numeric = numeric_distance,
      categorical = categorical_distance,
      interaction_active_average = active_average,
      interaction_level_delta = list(C1 = level_delta),
      interaction_observation = list(C1 = observation_delta)
    ),
    diagnostics = list(active_count = 1L, active_names = "C1"),
    gate = list(
      detected = c(C1 = TRUE),
      statistic = c(C1 = 0.7),
      p_value = c(C1 = 0.01)
    )
  )
  cache <- list(
    fits_by_prop = list(`0.1` = fit),
    default_gamma = 1,
    seconds_by_prop = c(`0.1` = 0)
  )
  setting <- data.frame(
    prop_nn = 0.1,
    gamma = 0.5,
    numeric_block_weight = 2
  )
  observed <- full_r1_distances_for_setting(cache, data.frame(x = 1:3), setting)
  expected_baseline <- 2 * numeric_distance + categorical_distance
  expected_added <- 0.5 * mean((2 * numeric_distance)[upper.tri(numeric_distance)]) *
    active_average

  testthat::expect_equal(
    observed$distances$Package_no_interaction_SC,
    expected_baseline,
    tolerance = 0
  )
  testthat::expect_equal(
    observed$distances$Gated_active_average_SC,
    expected_baseline + expected_added,
    tolerance = 1e-15
  )

  null_fit <- fit
  null_fit$components$interaction_active_average[,] <- 0
  null_fit$diagnostics$active_count <- 0L
  null_fit$diagnostics$active_names <- character()
  null_fit$gate$detected[["C1"]] <- FALSE
  cache$fits_by_prop[["0.1"]] <- null_fit
  null_observed <- full_r1_distances_for_setting(
    cache,
    data.frame(x = 1:3),
    setting
  )
  testthat::expect_identical(
    null_observed$distances$Gated_active_average_SC,
    null_observed$distances$Package_no_interaction_SC
  )
})

testthat::test_that("YAML parameter names and literal null state survive parsing", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  tasks <- full_expand_tasks(design)
  parameter_names <- unique(unlist(lapply(tasks, function(x) names(x$parameters))))

  testthat::expect_true("n" %in% parameter_names)
  testthat::expect_false(any(parameter_names %in% c("FALSE", "TRUE")))
  testthat::expect_false(any(grepl("^Var[0-9]+$", parameter_names)))

  states <- unique(unlist(lapply(tasks, function(x) x$parameters$interaction_state)))
  states <- states[!is.na(states)]
  testthat::expect_true(all(c("signal", "null") %in% states))
})

testthat::test_that("outer and within-task grids have the intended dimensions", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  tasks <- full_expand_tasks(design)
  counts <- table(vapply(tasks, `[[`, character(1), "scenario"))

  testthat::expect_equal(unname(counts[["corrected_marginal_overlap"]]), 150L)
  testthat::expect_equal(unname(counts[["categorical_contamination"]]), 200L)
  testthat::expect_equal(unname(counts[["pi_nn_sensitivity"]]), 100L)
  testthat::expect_equal(unname(counts[["gamma_sensitivity"]]), 100L)
  testthat::expect_equal(unname(counts[["numeric_block_weight_sensitivity"]]), 100L)
  testthat::expect_equal(unname(counts[["submitted_sim1_balanced"]]), 100L)
  testthat::expect_equal(unname(counts[["submitted_sim1_categorical_dominant"]]), 100L)
  testthat::expect_equal(unname(counts[["submitted_sim1_numeric_dominant"]]), 100L)
  testthat::expect_equal(unname(counts[["submitted_sim1_numeric_noise"]]), 100L)
  testthat::expect_equal(unname(counts[["submitted_sim2_interaction"]]), 100L)

  pi_task <- tasks[[which(vapply(tasks, `[[`, character(1), "scenario") ==
    "pi_nn_sensitivity")[1L]]]
  gamma_task <- tasks[[which(vapply(tasks, `[[`, character(1), "scenario") ==
    "gamma_sensitivity")[1L]]]
  weight_task <- tasks[[which(vapply(tasks, `[[`, character(1), "scenario") ==
    "numeric_block_weight_sensitivity")[1L]]]

  testthat::expect_equal(pi_task$within_task_grid$prop_nn, c(0.02, 0.05, 0.10, 0.20))
  testthat::expect_equal(
    gamma_task$within_task_grid$gamma,
    c(0, 0.25, 0.50, 0.75, 1)
  )
  testthat::expect_equal(
    weight_task$within_task_grid$numeric_block_weight,
    c(0.5, 1, 2)
  )
})

testthat::test_that("conditions are paired by scenario and replicate seed", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  tasks <- full_expand_tasks(design)
  overlap <- tasks[vapply(tasks, `[[`, character(1), "scenario") ==
    "corrected_marginal_overlap"]
  replicate_one <- overlap[vapply(overlap, `[[`, integer(1), "replicate") == 1L]

  testthat::expect_length(unique(vapply(replicate_one, `[[`, integer(1), "data_seed")), 1L)
  testthat::expect_equal(
    sort(unique(vapply(
      replicate_one,
      function(x) as.numeric(x$parameters$adjacent_overlap),
      numeric(1)
    ))),
    c(0.10, 0.25, 0.40)
  )
})

testthat::test_that("every manifest condition generates a valid smoke dataset", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  tasks <- full_expand_tasks(design)
  condition_key <- vapply(tasks, function(x) {
    paste(x$scenario, x$condition_id, sep = "::")
  }, character(1))
  representative_tasks <- tasks[!duplicated(condition_key)]

  for (task in representative_tasks) {
    generated <- full_generate_scenario(task, design, smoke = TRUE)
    testthat::expect_s3_class(generated$x, "data.frame")
    testthat::expect_equal(nrow(generated$x), length(generated$truth))
    testthat::expect_false(anyNA(generated$x))
    testthat::expect_true(all(vapply(
      generated$x,
      function(x) is.numeric(x) || is.factor(x),
      logical(1)
    )))
    testthat::expect_gte(nlevels(factor(generated$truth)), 2L)
  }
})

testthat::test_that("submitted-simulation continuity generators preserve their declared designs", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  tasks <- full_expand_tasks(design)
  first_task <- function(scenario) {
    tasks[[which(vapply(tasks, `[[`, character(1), "scenario") == scenario)[1L]]]
  }

  simulation_1 <- full_generate_scenario(
    first_task("submitted_sim1_balanced"),
    design,
    smoke = TRUE
  )
  testthat::expect_equal(dim(simulation_1$x), c(90L, 30L))
  testthat::expect_equal(as.integer(table(simulation_1$truth)), c(18L, 27L, 45L))
  testthat::expect_true(all(vapply(simulation_1$x[16:30], is.factor, logical(1))))
  testthat::expect_gt(
    simulation_1$diagnostics$covariance[[3L]]$clipped_eigenvalues,
    0L
  )

  simulation_2 <- full_generate_scenario(
    first_task("submitted_sim2_interaction"),
    design,
    smoke = TRUE
  )
  testthat::expect_equal(dim(simulation_2$x), c(90L, 9L))
  testthat::expect_true(all(vapply(simulation_2$x[7:9], is.factor, logical(1))))
  testthat::expect_equal(as.integer(table(simulation_2$x$C1)), c(36L, 18L, 36L))
  testthat::expect_equal(as.integer(table(simulation_2$x$C2)), c(72L, 18L))
  testthat::expect_equal(as.integer(table(simulation_2$x$C3)), c(54L, 36L))
})

testthat::test_that("legacy Simulation 1 generator preserves the executable design", {
  testthat::skip_if_not_installed("mvtnorm")
  parameters <- list(n = 100L, p_numeric = 15L, q_categorical = 15L)

  observed <- full_make_legacy_simulation_1(parameters, seed = 731L)
  repeated <- full_make_legacy_simulation_1(parameters, seed = 731L)

  testthat::expect_identical(observed, repeated)
  testthat::expect_equal(dim(observed$x), c(100L, 30L))
  testthat::expect_equal(as.integer(table(observed$truth)), c(20L, 30L, 50L))
  testthat::expect_true(all(vapply(observed$x[16:30], is.factor, logical(1))))
  testthat::expect_equal(
    unname(vapply(observed$x[16:30], nlevels, integer(1))),
    rep(2L, 15L)
  )
  testthat::expect_identical(
    observed$diagnostics$source,
    "historical_executable_code"
  )
  testthat::expect_true(
    observed$diagnostics$categorical_generated_before_numeric_shifts
  )
  testthat::expect_gt(
    observed$diagnostics$covariance[[3L]]$clipped_eigenvalues,
    0L
  )
})

testthat::test_that("legacy Simulation 2 generator preserves counts and recursion", {
  parameters <- list(n = 100L, p_numeric = 6L, q_categorical = 4L)
  observed <- full_make_legacy_simulation_2(parameters, seed = 913L)
  repeated <- full_make_legacy_simulation_2(parameters, seed = 913L)

  testthat::expect_identical(observed, repeated)
  testthat::expect_equal(dim(observed$x), c(100L, 10L))
  testthat::expect_equal(as.integer(table(observed$x$C1)), c(40L, 20L, 40L))
  testthat::expect_equal(as.integer(table(observed$x$C2)), c(80L, 20L))
  testthat::expect_equal(as.integer(table(observed$x$C3)), c(60L, 40L))
  testthat::expect_equal(as.integer(table(observed$x$C4)), c(40L, 30L, 30L))

  coefficient_matrix <- observed$diagnostics$coefficient_matrix
  coefficient <- coefficient_matrix[cbind(
    as.integer(observed$x$C1),
    as.integer(observed$x$C2)
  )]
  intercept <- ifelse(coefficient == 0.5, 8, 0)
  for (variable in 3:1) {
    testthat::expect_equal(
      observed$x[[paste0("V", variable)]],
      intercept + coefficient * rowSums(
        observed$x[paste0("V", (variable + 1L):6L)]
      ),
      tolerance = 1e-12
    )
  }
  expected_truth <- factor(
    ifelse(coefficient == 0.8, "C1", ifelse(coefficient == -0.4, "C2", "C3")),
    levels = c("C1", "C2", "C3")
  )
  testthat::expect_identical(observed$truth, expected_truth)
})

testthat::test_that("legacy Simulation 2 rejects dimensions absent from the source code", {
  testthat::expect_error(
    full_make_legacy_simulation_2(
      list(n = 100L, p_numeric = 6L, q_categorical = 3L),
      seed = 1L
    ),
    "fixes `p_numeric = 6` and `q_categorical = 4`",
    fixed = TRUE
  )
})

testthat::test_that("irrelevant-variable grids preserve the fixed signal draws", {
  design <- full_read_design(file.path(full_dir, "design.yml"))
  tasks <- full_expand_tasks(design)

  first_replicate <- function(scenario) {
    tasks[
      vapply(tasks, `[[`, character(1), "scenario") == scenario &
        vapply(tasks, `[[`, integer(1), "replicate") == 1L
    ]
  }

  numeric_tasks <- first_replicate("irrelevant_numeric")
  numeric_data <- lapply(
    numeric_tasks,
    full_generate_scenario,
    design = design,
    smoke = FALSE
  )
  numeric_reference <- numeric_data[[1L]]
  for (generated in numeric_data[-1L]) {
    testthat::expect_identical(generated$truth, numeric_reference$truth)
    testthat::expect_identical(
      generated$x[grep("^V_signal_", names(generated$x))],
      numeric_reference$x[grep("^V_signal_", names(numeric_reference$x))]
    )
    testthat::expect_identical(
      generated$x[grep("^C_signal_", names(generated$x))],
      numeric_reference$x[grep("^C_signal_", names(numeric_reference$x))]
    )
  }

  categorical_tasks <- first_replicate("irrelevant_categorical")
  categorical_data <- lapply(
    categorical_tasks,
    full_generate_scenario,
    design = design,
    smoke = FALSE
  )
  categorical_reference <- categorical_data[[1L]]
  for (generated in categorical_data[-1L]) {
    testthat::expect_identical(generated$truth, categorical_reference$truth)
    testthat::expect_identical(
      generated$x[grep("^V_signal_", names(generated$x))],
      categorical_reference$x[grep("^V_signal_", names(categorical_reference$x))]
    )
    testthat::expect_identical(
      generated$x[grep("^C_signal_", names(generated$x))],
      categorical_reference$x[grep("^C_signal_", names(categorical_reference$x))]
    )
  }
})
