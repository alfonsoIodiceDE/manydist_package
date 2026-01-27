# data-raw/mdist_methods_info.R
mdist_param_sets <- list(
  baseline = list(
    preset        = "custom",
    distance_cont = "manhattan",
    commensurable = FALSE,
    scaling_cont  = "std",
    ncomp         = NULL,
    threshold     = NULL
  ),
  hl = list(
    preset        = "custom",
    distance_cont = "euclidean",
    distance_cat  = "HLeucl",
    commensurable = FALSE,
    scaling_cont  = "std",
    ncomp         = NULL,
    threshold     = NULL
  ),
  hl_add = list(
    preset        = "custom",
    distance_cont = "manhattan",
    distance_cat  = "HL",
    commensurable = FALSE,
    scaling_cont  = "std",
    ncomp         = NULL,
    threshold     = NULL
  ),
  u_ind = list(
    preset        = "custom",
    distance_cont = "manhattan",
    distance_cat  = "matching",
    commensurable = TRUE,
    scaling_cont  = "std",
    ncomp         = NULL,
    threshold     = NULL
  ),
  u_mix = list(
    preset        = "custom",
    distance_cont = "manhattan",
    distance_cat  = "tot_var_dist",
    commensurable = TRUE,
    scaling_cont  = "std",
    ncomp         = NULL,
    threshold     = NULL
  ),
  u_dep_eucl = list(
    preset        = "custom",
    distance_cont = "euclidean",
    distance_cat  = "tot_var_dist",
    commensurable = TRUE,
    scaling_cont  = "std",
    ncomp         = NULL,
    threshold     = NULL
  )
)

mdist_method_lookup <- tibble::tibble(
  method = c(
    "baseline", "naive", "gower", "hl", "hl_add",
    "u_dep", "u_ind", "u_mix", "u_dep_eucl",
    "gudmm", "dkss", "mod_gower"
    # consider dropping "euclidean_onehot" or "naive" to avoid duplication
  ),
  mdist_type = c(
    "custom", "preset", "preset", "custom", "custom",
    "preset", "custom", "custom", "custom",
    "preset", "preset", "preset"
  ),
  mdist_preset = c(
    NA_character_, "euclidean_onehot", "gower", NA_character_, NA_character_,
    "unbiased_dependent", NA_character_, NA_character_, NA_character_,
    "gudmm", "dkss", "mod_gower"
  ),
  param_set = list(
    mdist_param_sets$baseline,
    NULL,
    NULL,
    mdist_param_sets$hl,
    mdist_param_sets$hl_add,
    NULL,
    mdist_param_sets$u_ind,
    mdist_param_sets$u_mix,
    mdist_param_sets$u_dep_eucl,
    NULL,
    NULL,
    NULL
  )
)

usethis::use_data(
  mdist_param_sets,
  mdist_method_lookup,
  internal = TRUE,
  overwrite = TRUE
)
