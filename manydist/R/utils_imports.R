#' @keywords internal
#' @noRd

# magrittr
#' @importFrom magrittr %>%
NULL

# dplyr / tidyselect / tidyr / tibble / purrr
#' @importFrom dplyr mutate rename select across pull
#' @importFrom tidyselect where everything matches any_of all_of
#' @importFrom tidyr drop_na crossing
#' @importFrom tibble tibble as_tibble is_tibble
#' @importFrom purrr map map2 map_dbl map2_dbl pmap map_df
NULL

# recipes / rsample
#' @importFrom recipes recipe prep bake step_dummy step_range step_center step_select
#' @importFrom recipes step_scale step_pca step_mutate all_predictors all_nominal all_nominal_predictors all_numeric_predictors
#' @importFrom rsample initial_split training testing
NULL

# cluster / stats / distances / Matrix
#' @importFrom cluster daisy
#' @importFrom stats dist as.dist cov var sd median density model.matrix IQR quantile filter start
#' @importFrom distances distances
#' @importFrom Matrix bdiag
NULL

# data.table / fastDummies / forcats / readr / entropy
#' @importFrom data.table as.data.table data.table setDT setnames rbindlist
#' @importFrom data.table `:=`
#' @importFrom fastDummies dummy_cols
#' @importFrom forcats fct_drop fct_count fct_lump_lowfreq
#' @importFrom readr parse_number
#' @importFrom entropy mi.empirical
NULL

# ggplot2 generic if you provide autoplot.MDistLOVO
#' @importFrom ggplot2 autoplot
NULL
