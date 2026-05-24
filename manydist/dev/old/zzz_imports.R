#' Package imports
#'
#' @name manydist-imports
#' @keywords internal
#' @noRd
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate rename select across pull filter arrange group_by summarise slice_min slice_max desc
#' @importFrom tidyselect where everything matches any_of all_of starts_with
#' @importFrom tidyr drop_na crossing
#' @importFrom tibble tibble as_tibble is_tibble
#' @importFrom purrr map map2 map_dbl map2_dbl pmap map_df map_chr
#' @importFrom recipes recipe prep bake step_dummy step_range step_center step_select step_unknown step_novel
#' @importFrom recipes step_scale step_pca step_mutate step_normalize step_rm
#' @importFrom recipes all_predictors all_nominal all_nominal_predictors all_numeric_predictors
#' @importFrom rsample initial_split training testing
#' @importFrom cluster daisy
#' @importFrom stats dist as.dist cov var sd median density model.matrix IQR quantile  start
#' @importFrom stats kmeans predict rchisq rnorm runif
#' @importFrom utils capture.output
#' @importFrom distances distances
#' @importFrom Matrix bdiag
#' @importFrom data.table as.data.table data.table setDT setnames rbindlist
#' @importFrom data.table `:=`
#' @importFrom fastDummies dummy_cols
#' @importFrom forcats fct_drop fct_count fct_lump_lowfreq
#' @importFrom readr parse_number
#' @importFrom entropy mi.empirical
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#' @importFrom tune tunable
NULL

utils::globalVariables(c(
  "Species",
  "argument",
  "data_type",
  "distance_basis",
  "engine",
  "method",
  "new_data",
  "object",
  "spec_type"
))
