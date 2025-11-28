#' Build a recipe that computes mdist-based distance features
#'
#' @description
#' This helper creates a tidymodels recipe using step_mdist().
#' It supports both mdist presets and custom param sets.
#'
#' @param df A data frame.
#' @param mdist_type "preset" or "custom".
#' @param mdist_preset Name of the preset (if mdist_type == "preset").
#' @param param_set A list of custom mdist arguments (if mdist_type == "custom").
#' @param outcome Name of the outcome variable.
#'
#' @return A `recipes::recipe()` object.
#'
#' @export
make_mdist_recipe <- function(df, mdist_type, mdist_preset, param_set, outcome) {
  # 1. Outcome must exist
  if (!outcome %in% names(df)) {
    rlang::abort(paste0("Outcome '", outcome, "' not found in data."))
  }

  # 2. There must be at least one predictor
  pred_names <- setdiff(names(df), outcome)
  if (length(pred_names) == 0L) {
    rlang::abort("No predictors: data only contains the outcome column.")
  }

  # 3. Base recipe: outcome ~ .
  rec <- recipes::recipe(
    stats::reformulate(termlabels = ".", response = outcome),
    data = df
  )

  # 4. Attach step_mdist depending on mdist_type -------------------------
  if (identical(mdist_type, "preset")) {

    rec <- rec |>
      step_mdist(
        all_of(pred_names),
        preset = mdist_preset
      )

  } else if (identical(mdist_type, "custom")) {

    if (is.null(param_set)) {
      rlang::abort("`mdist_type = 'custom'` but `param_set` is NULL.")
    }

    # param_set must be a *named list* with entries like:
    # list(preset="custom", distance_cont="euclidean", ...)
    args <- c(
      list(               # first two arguments to step_mdist()
        recipe = rec,
        all_of(pred_names)
      ),
      param_set           # named args: preset, distance_cont, ...
    )

    rec <- do.call(step_mdist, args)

  } else if (identical(mdist_type, "none") || is.na(mdist_type)) {
    # no mdist step: leave predictors as is
    rec <- rec

  } else {
    rlang::abort(paste0("Unknown mdist_type: ", mdist_type))
  }

  rec
}
