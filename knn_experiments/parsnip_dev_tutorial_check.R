library(tidymodels)

## model + modes ---------------------------------------------------------
set_new_model("nearest_neighbor_dist")

set_model_mode("nearest_neighbor_dist", "classification")
set_model_mode("nearest_neighbor_dist", "regression")

## engine name: choose "precomputed" for distances
set_model_engine(
  "nearest_neighbor_dist",
  mode = "classification",
  eng  = "precomputed"
)

set_model_engine(
  "nearest_neighbor_dist",
  mode = "regression",
  eng  = "precomputed"
)

## if this is in the package, you'd also do:
 set_dependency("nearest_neighbor_dist", eng = "precomputed", pkg = "manydist")

## arguments -------------------------------------------------------------
set_model_arg(
  model    = "nearest_neighbor_dist",
  eng      = "precomputed",
  parsnip  = "neighbors",   # parsnip argument name
  original = "k",           # argument name in fit_knn_dist()
  func     = list(pkg = "rlang", fun = "quo"),
  has_submodel = FALSE
)

set_model_arg(
  model    = "nearest_neighbor_dist",
  eng      = "precomputed",
  parsnip  = "dist_fun",
  original = "dist_fun",
  func     = list(pkg = "rlang", fun = "quo"),
  has_submodel = FALSE
)

set_model_arg(
  model    = "nearest_neighbor_dist",
  eng      = "precomputed",
  parsnip  = "dist_args",
  original = "dist_args",
  func     = list(pkg = "rlang", fun = "quo"),
  has_submodel = FALSE
)

show_model_info("nearest_neighbor_dist")


nearest_neighbor_dist <-
  function(mode      = "classification",
           neighbors = NULL,
           dist_fun  = NULL,
           dist_args = NULL) {
    if (!mode %in% c("classification", "regression", "unknown")) {
      rlang::abort("`mode` must be 'classification', 'regression', or 'unknown'.")
    }

    args <- list(
      neighbors = rlang::enquo(neighbors),
      dist_fun  = rlang::enquo(dist_fun),
      dist_args = rlang::enquo(dist_args)
    )

    new_model_spec(
      "nearest_neighbor_dist",
      args     = args,
      eng_args = NULL,
      mode     = mode,
      method   = NULL,
      engine   = NULL
    )
  }


set_fit(
  model = "nearest_neighbor_dist",
  eng   = "precomputed",
  mode  = "classification",
  value = list(
    interface = "data.frame",
    protect   = c("x", "y"),
    func      = c(pkg = NULL, fun = "fit_knn_dist"),
    defaults  = list()
  )
)

set_fit(
  model = "nearest_neighbor_dist",
  eng   = "precomputed",
  mode  = "regression",
  value = list(
    interface = "data.frame",
    protect   = c("x", "y"),
    func      = c(pkg = NULL, fun = "fit_knn_dist"),
    defaults  = list()
  )
)

show_model_info("nearest_neighbor_dist")

nearest_neighbor_dist(neighbors = 5, dist_fun = dummy_dist) |>
  translate(engine = "precomputed")



dat <- as_tibble(iris) |>
  mutate(Species = factor(Species))

set.seed(1)
sp  <- initial_split(dat, prop = 0.7, strata = Species)
tr  <- training(sp)
te  <- testing(sp)

spec <- nearest_neighbor_dist(
  mode      = "classification",
  neighbors = 5,
  dist_fun  = dummy_dist,
  dist_args = list()          # you can leave this out, has default in fit_knn_dist
) |>
  set_engine("precomputed")


fit_parsnip <- fit(
  spec,
  Species ~ .,              # formula: y ~ predictors
  data = tr
)

str(fit_parsnip$fit)        # should look like your knn_dist object
