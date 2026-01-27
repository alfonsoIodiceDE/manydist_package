rm(list=ls())
packages <- c(
  "aricode", "cluster", "fpc", "manydist", "patchwork", "tidyverse",
  "tidymodels", "varhandle", "vegan", "kdml", "ggplot2", "viridis",
  "dplyr", "entropy", "kernlab", "mclust", "Matrix", "philentropy", "conflicted"
)

# Install any missing packages from CRAN
missing <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(missing) > 0) {
  install.packages(missing)
}

# Load all packages quietly
invisible(
  lapply(packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  })
)

# Optional: Resolve common function conflicts
# install.packages("conflicted")  # Only once
tidymodels_prefer()
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("distance", "philentropy")

reps<-100  # 100 can take 30min
n<-500
k<-2 # Dimension solution
sigma<-0.03 # Noise

porig<-6 # Number of variables corresponding to true config
pnnoise<-0 # Number of numerical noise vars
pcnoise<-0 # Number of categorical noise vars
pn<-2 # Number of numerical variables underying config

pnum<-pn+pnnoise+pcnoise  # Total number of numerical vars+noise vars (Needed to know which vars to discretize)
p<-porig+pnnoise+pcnoise # Total number of variables
qoptions<-c(2,3,5,9) # Distribution of the categorical variables corresponding to true config
pcat<-length(qoptions) # Number of cat variables underlying config


methods = c("baseline","naive","hl","hl_add","gower","gudmm","dkss","mod_gower","u_ind","u_dep","u_mix")
method_levels <- c("baseline", "naive", "hl", "hl_add", "gower","gudmm", "dkps", "mod_gower", "u_ind", "u_dep", "u_mix")


my_colors <- c(
  # Set B
  "baseline"    = "#8C564B",  # brown
  "naive"       = "#E377C2",  # pink
  "hl"          = "#9C9E00",  # darker yellow-green (better contrast)
  "gudmm"       = "#FF7F0E",  # orange
  "dkps"        = "#7F7F7F",  # darker grey (more visible)

  # Set A
  "hl_add"      = "#9EC5E5",  # light blue
  "u_ind"       = "#2CA02C",  # green
  "u_dep"       = "#D62728",  # red
  "u_mix"       = "#0D3B66",  # dark blue
  "gower"       = "#FFBF00",  # gold
  "mod_gower"   = "#9467BD"   # purple
)

my_shapes <- c(
  "baseline"    = 16,  # filled circle
  "naive"       = 17,  # filled triangle
  "hl"          = 15,  # filled square
  "hl_add"      = 3,   # plus
  "gower"       = 4,   # cross
  "gudmm"       = 18,  # filled diamond
  "dkps"        = 8,   # star
  "mod_gower"   = 0,   # empty square
  "u_ind"       = 1,   # empty circle
  "u_dep"       = 2,    # empty triangle
  "u_mix" = 5 # empty diamond
)

library(progress)


pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent eta: :eta",
  total = reps * length(methods),
  clear = FALSE,
  width = 60
)

# methods = c("baseline","naive","gower","hl","u_dep","u_ind","gudmm","dkss","mod_gower")


generated_datasets = tibble(replicate=1:reps) |>
  mutate(dataset = map(.x=replicate, ~generate_dataset(n = n, pn=pn, porig=porig, pnnoise =pnnoise , pcnoise = pcnoise, sigma = sigma, qoptions = qoptions, seed = 1234+.x)),
         method=map(.x=replicate,~methods)
  ) |> unnest(method)|>
  left_join(mdist_method_lookup, by = "method")



run_mdist_within <- function(df, mdist_type, mdist_preset, param_set, outcome = "response") {
  x <- df

  if (is.character(outcome) && length(outcome) == 1 && outcome %in% names(df)) {
    x <- dplyr::select(df, -dplyr::all_of(outcome))
  }

  if (mdist_type == "preset") {
    out <- manydist::mdist(x = x, preset = mdist_preset)

  } else if (mdist_type == "custom") {
    if (is.null(param_set)) stop("param_set is NULL for mdist_type = 'custom'")
    args <- param_set
    args$x <- x
    out <- do.call(manydist::mdist, args = args)

  } else {
    stop("Unknown mdist_type: ", mdist_type)
  }

  out$distance
}

# little helper for the leave one variable out version of mdist
run_lovo <- function(df, mdist_type, mdist_preset, param_set,
                     outcome = "response", dims = 2, keep_dist = FALSE) {

  x <- df  # default: don't drop anything

  if (is.character(outcome) && length(outcome) == 1 && outcome %in% names(df)) {
    x <- dplyr::select(df, -dplyr::all_of(outcome))
  }

  if (mdist_type == "preset") {
    manydist::lovo_mdist(x = x, preset = mdist_preset, dims = dims, keep_dist = keep_dist)

  } else if (mdist_type == "custom") {
    args <- param_set
    args$x <- x
    args$dims <- dims
    args$keep_dist <- keep_dist
    do.call(manydist::lovo_mdist, args = args)

  } else {
    stop("Unknown mdist_type: ", mdist_type)
  }
}


simulation_structure <- generated_datasets |>
  dplyr::mutate(
    loo_results = purrr::pmap(
      dplyr::pick(dataset, method, mdist_type, mdist_preset, param_set),
      \(dataset, method, mdist_type, mdist_preset, param_set) {
        pb$tick()
        X <- if (method == "baseline") dataset$Xorig else dataset$X
        run_lovo(X, mdist_type, mdist_preset, param_set)
      }
    )
  )


save(file="simulation_structure.RData", simulation_structure)

recovery_structure <- tidyr::crossing(
  tibble::tibble(replicate = 1:reps),
  method = methods,
  categories = qoptions
) |>
  dplyr::left_join(manydist:::mdist_method_lookup, by = "method") |>
  dplyr::mutate(
    dataset = purrr::map2(
      .x = replicate, .y = categories,
      ~ generate_dataset(
        n = n, porig = porig, pn = pn,
        pnnoise = pnnoise, pcnoise = pcnoise,
        sigma = sigma, qoptions = .y, mode = "shared",
        seed = 1234 + .x
      )
    ),
    # compute distance matrix using method specs + baseline uses Xorig
    distance_matrix = purrr::pmap(
      dplyr::pick(dataset, method, mdist_type, mdist_preset, param_set),
      \(dataset, method, mdist_type, mdist_preset, param_set) {
        X <- if (method == "baseline") dataset$Xorig else dataset$X

        run_mdist_within(
          df = X,
          mdist_type = mdist_type,
          mdist_preset = mdist_preset,
          param_set = param_set,
          outcome = NULL   # because X / Xorig have no response column
        ) |>
          base::as.matrix()
      }
    ),

    mds_results = purrr::map(
      distance_matrix,
      ~ cmdscale(.x, eig = TRUE, k = 2)$points |>
        as.data.frame() |>
        setNames(c("x1", "x2"))
    ),

    congruence_coeff = purrr::map2(
      .x = mds_results, .y = dataset,
      \(mds, ds) congruence_coeff(ds$truth, mds)
    )
  )

# save(file="recovery_structure.RData", recovery_structure)
recovery_structure_lite = recovery_structure |>
  select(-dataset, -distance_matrix, -mds_results) |>
  unnest(congruence_coeff) |>
  select(replicate, method, categories, congruence_coeff) |>
  mutate(categories=factor(paste0(categories," categories")),
         alienation_coeff = sqrt(1-congruence_coeff^2))

save(file="recovery_structure_lite.RData", recovery_structure_lite)



k_true <- 4
q <- 9
numsep <- 0.1
catsep <- 0.6
clustSizeEq <- 50

param_generator_grid <- tibble::tibble(
  numsignal = c(0,0,4,4,4,8,8,8),
  numnoise  = c(8,8,4,4,4,0,0,0),
  catsignal = c(4,8,0,4,8,0,4,8),
  catnoise  = c(4,0,8,4,0,8,4,0)
) |>
  tidyr::crossing(replicate = 1:reps)

mixed_datasets_grid <- param_generator_grid |>
  dplyr::mutate(
    data = purrr::pmap(
      dplyr::pick(numsignal, catsignal, numnoise, catnoise, replicate),
      \(numsignal, catsignal, numnoise, catnoise, replicate) {
        gen_mixed(
          k_true      = k_true,
          clustSizeEq = clustSizeEq,
          numsignal   = numsignal,
          numnoise    = numnoise,
          catsignal   = catsignal,
          catnoise    = catnoise,
          q           = q,
          q_err       = q,
          numsep      = numsep,
          catsep      = catsep,
          seed        = replicate
        )
      }
    )
  )

methods <- c("naive","hl","hl_add","gower","gudmm","dkss","mod_gower","u_ind","u_dep","u_mix")

sim_methods_grid <- mixed_datasets_grid |>
  tidyr::crossing(method = methods) |>
  dplyr::left_join(manydist:::mdist_method_lookup, by = "method")

stopifnot(!any(is.na(sim_methods_grid$mdist_type)))

ari_pam_results <- sim_methods_grid |>
  dplyr::mutate(
    ari = purrr::pmap_dbl(
      dplyr::pick(data, mdist_type, mdist_preset, param_set),
      \(data, mdist_type, mdist_preset, param_set) {
        df <- data$df
        D <- run_mdist_within(df, mdist_type, mdist_preset, param_set, outcome = "y")
        pam_fit <- cluster::pam(D, k = k_true, diss = TRUE)
        mclust::adjustedRandIndex(df$y, pam_fit$clustering)
      }
    )
  )

save(ari_pam_results, file = "ari_pam_experiment.RData")
