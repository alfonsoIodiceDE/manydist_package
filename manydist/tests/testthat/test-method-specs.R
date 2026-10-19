test_that("preset selection replaces mode without including unwanted rows", {
  expect_false("mode" %in% names(formals(all_dist_method_specs)))
  expect_error(all_dist_method_specs(mode = "full"), "unused argument")
  custom <- all_dist_method_specs()
  expect_true(all(custom$preset == "custom"))
  expect_false(anyNA(custom))

  selected <- all_dist_method_specs(preset = c("gower", "u_indep", "u_dep"))
  expect_setequal(selected$preset, c("gower", "u_indep", "u_dep"))
  expect_equal(nrow(selected), 3L)
  expect_true(all(is.na(selected$method_cat)))

  presets <- all_dist_method_specs(preset = "all")
  catalogue <- dist_methods_tbl() |>
    dplyr::filter(argument == "preset", method != "custom")
  expect_setequal(presets$preset, catalogue$method)
  expect_false("custom" %in% presets$preset)
  expect_identical(all_dist_method_specs(preset = "custom"), custom)
  full <- all_dist_method_specs(preset = c("all", "custom"))
  expect_equal(nrow(full), nrow(presets) + nrow(custom))
  expect_identical(dplyr::filter(full, preset == "custom"), custom)
  expect_equal(nrow(all_dist_method_specs(preset = c("gower", "gower"))), 1L)
  expect_equal(nrow(all_dist_method_specs(preset = character())), 0L)
})

test_that("component arguments filter custom grids and can accompany named presets", {
  selected <- all_dist_method_specs(
    method_cat = c("matching", "tvd"),
    method_num = c("std", "robust", "pc_scores"), commensurable = TRUE
  )
  expect_equal(nrow(selected), 4L)
  expect_setequal(selected$method_num, c("std", "pc_scores"))
  expect_true(all(selected$commensurable))

  mixed <- all_dist_method_specs(
    preset = c("gower", "custom"), method_cat = "matching",
    method_num = "std", commensurable = TRUE
  )
  expect_equal(nrow(mixed), 2L)
  expect_setequal(mixed$preset, c("gower", "custom"))
  expect_true(is.na(mixed$commensurable[mixed$preset == "gower"]))

  unscaled <- all_dist_method_specs(method_cat = "matching", commensurable = FALSE)
  expect_setequal(unscaled$method_num, c("none", "std", "robust", "range", "pc_scores"))
  expect_equal(nrow(all_dist_method_specs(method_cat = character())), 0L)
  expect_equal(nrow(all_dist_method_specs(commensurable = logical())), 0L)
})

test_that("specification selection rejects invalid or silently ignored filters", {
  expect_error(all_dist_method_specs(preset = "unknown"), "`preset`")
  expect_error(all_dist_method_specs(preset = NA_character_), "`preset`")
  expect_error(all_dist_method_specs(preset = TRUE), "`preset`")
  expect_error(all_dist_method_specs(method_num = "unknown"), "`method_num`")
  expect_error(all_dist_method_specs(method_cat = NA_character_), "`method_cat`")
  expect_error(all_dist_method_specs(commensurable = NA), "`commensurable`")
  expect_error(all_dist_method_specs(commensurable = 1), "`commensurable`")
  expect_error(all_dist_method_specs(preset = "all", commensurable = TRUE),
               "Component filters apply to custom rows")
})

test_that("affine representatives remain equivalent with fitted training means", {
  train <- data.frame(z = c(0, 2, 8, 10), w = c(1, 5, 2, 7))
  test <- data.frame(z = c(1, 100), w = c(3, -100))
  for (new_data in list(NULL, test)) {
    standard <- mdist(train, new_data = new_data, method_num = "std")
    robust <- mdist(train, new_data = new_data, method_num = "robust")
    expect_equal(standard$distance, robust$distance)
  }
})
