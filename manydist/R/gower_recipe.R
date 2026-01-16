gower_recipe <- function(data) {
  .x <- NULL
  n_cat <- data |> dplyr::select(where(is.factor)) |> ncol()

  if (n_cat > 0) {
    gow_recipe <- recipe(~., data = data) |>
      step_range(all_numeric_predictors(), min = 0, max = 1) |>
      step_unknown(all_nominal_predictors()) |>
      step_novel(all_nominal_predictors()) |>
      step_dummy(
        all_nominal_predictors(),
        one_hot = TRUE, id = "dummy",
        naming = function(var, lvl, ordinal, ...) paste0("dummy_", var, "_", lvl)
      ) |>
      step_mutate(across(starts_with("dummy_"), ~ .x / 2))
  } else {
    gow_recipe <- recipe(~., data = data) |>
      step_range(all_numeric_predictors(), min = 0, max = 1)
  }

  gow_recipe
}
