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
  u_dep_manh = list(
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

mdist_method_lookup <- tibble(
  method = c(
    "baseline",
    "naive",
    "gower",
    "hl",
    "hl_add",
    "u_dep",
    "u_ind",
    "u_dep_manh",
    "u_dep_eucl",
    "gudmm",
    "dkss",
    "mod_gower",
    "euclidean_onehot"  # if you use method name = preset name
  ),
  mdist_type = c(
    "custom",    # baseline
    "preset",    # naive -> euclidean_onehot preset
    "preset",    # gower
    "custom",    # hl
    "custom",    # hl_add
    "preset",    # u_dep -> unbiased_dependent preset
    "custom",    # u_ind
    "custom",    # u_dep_manh
    "custom",    # u_dep_eucl
    "preset",    # gudmm
    "preset",    # dkss
    "preset",    # mod_gower
    "preset"     # euclidean_onehot
  ),
  mdist_preset = c(
    NA_character_,          # baseline
    "euclidean_onehot",     # naive
    "gower",                # gower
    NA_character_,          # hl
    NA_character_,          # hl_add
    "unbiased_dependent",   # u_dep
    NA_character_,          # u_ind
    NA_character_,          # u_dep_manh
    NA_character_,          # u_dep_eucl
    "gudmm",                # gudmm
    "dkss",                 # dkss
    "mod_gower",            # mod_gower
    "euclidean_onehot"      # euclidean_onehot
  ),
  param_set = list(
    mdist_param_sets$baseline,  # baseline
    NULL,                       # naive
    NULL,                       # gower
    mdist_param_sets$hl,        # hl
    mdist_param_sets$hl_add,    # hl_add
    NULL,                       # u_dep (preset)
    mdist_param_sets$u_ind,     # u_ind
    mdist_param_sets$u_dep_manh,# u_dep_manh
    mdist_param_sets$u_dep_eucl,# u_dep_eucl
    NULL,                       # gudmm
    NULL,                       # dkss
    NULL,                       # mod_gower
    NULL                        # euclidean_onehot
  )
)

