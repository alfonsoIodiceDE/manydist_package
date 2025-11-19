## =====================================================================
##  TEST SCRIPT: nearest_neighbor_dist + mdist + tidymodels
## =====================================================================

# Clean up ----------------------------------------------------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(tidymodels)
  library(pkgload)
  library(dplyr)
  library(tibble)
})

## -----------------------------------------------------------------------
## 1. Load development version of manydist
## -----------------------------------------------------------------------

# Adjust this path if your script lives somewhere else
pkg_path <- "../manydist"

if ("manydist" %in% loadedNamespaces()) {
  unloadNamespace("manydist")
}

pkgload::load_all(pkg_path, export_all = TRUE, helpers = TRUE, attach_testthat = FALSE)

cat("Loaded manydist from:\n  ",
    getNamespaceInfo(asNamespace("manydist"), "path"), "\n\n", sep = "")


# ------------------------------------------------------------------
# DEV-ONLY registration: use pkg = NULL so functions are found in
# the global environment (via export_all=TRUE)
# ------------------------------------------------------------------
register_nearest_neighbor_dist_model_dev <- function() {
  q <- function(expr) try(suppressWarnings(expr), silent = TRUE)

  # model + modes
  q(parsnip::set_new_model("nearest_neighbor_dist"))
  q(parsnip::set_model_mode("nearest_neighbor_dist", "classification"))
  q(parsnip::set_model_mode("nearest_neighbor_dist", "regression"))

  # engine
  q(parsnip::set_model_engine("nearest_neighbor_dist",
                              mode = "classification",
                              eng  = "precomputed"))
  q(parsnip::set_model_engine("nearest_neighbor_dist",
                              mode = "regression",
                              eng  = "precomputed"))

  # args
  q(parsnip::set_model_arg(
    model    = "nearest_neighbor_dist",
    eng      = "precomputed",
    parsnip  = "neighbors",
    original = "neighbors",
    func     = list(pkg = "rlang", fun = "quo"),
    has_submodel = FALSE
  ))
  q(parsnip::set_model_arg(
    model    = "nearest_neighbor_dist",
    eng      = "precomputed",
    parsnip  = "dist_fun",
    original = "dist_fun",
    func     = list(pkg = "rlang", fun = "quo"),
    has_submodel = FALSE
  ))
  q(parsnip::set_model_arg(
    model    = "nearest_neighbor_dist",
    eng      = "precomputed",
    parsnip  = "dist_args",
    original = "dist_args",
    func     = list(pkg = "rlang", fun = "quo"),
    has_submodel = FALSE
  ))

  # FIT — note pkg = NULL here
  q(parsnip::set_fit(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "classification",
    value = list(
      interface = "data.frame",
      protect   = c("x", "y"),
      func      = c(pkg = NULL, fun = "fit_nearest_neighbor_dist"),
      defaults  = list()
    )
  ))
  q(parsnip::set_fit(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "regression",
    value = list(
      interface = "data.frame",
      protect   = c("x", "y"),
      func      = c(pkg = NULL, fun = "fit_nearest_neighbor_dist"),
      defaults  = list()
    )
  ))

  # PREDICT — note pkg = NULL here too
  q(parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "classification",
    type  = "class",
    value = list(
      pre  = NULL, post = NULL,
      func = c(pkg = NULL, fun = "predict_nearest_neighbor_dist_class"),
      args = list(
        object   = quote(object$fit),
        new_data = quote(new_data),
        type     = "class"
      )
    )
  ))
  q(parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "classification",
    type  = "prob",
    value = list(
      pre  = NULL, post = NULL,
      func = c(pkg = NULL, fun = "predict_nearest_neighbor_dist_class"),
      args = list(
        object   = quote(object$fit),
        new_data = quote(new_data),
        type     = "prob"
      )
    )
  ))
  q(parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL, post = NULL,
      func = c(pkg = NULL, fun = "predict_nearest_neighbor_dist_reg"),
      args = list(
        object   = quote(object$fit),
        new_data = quote(new_data)
      )
    )
  ))

  invisible(TRUE)
}

# Call the DEV registration
register_nearest_neighbor_dist_model_dev()
parsnip::show_model_info("nearest_neighbor_dist")

## -----------------------------------------------------------------------
## 2. Register the parsnip model and inspect it
## -----------------------------------------------------------------------

# instead of manydist::register_nearest_neighbor_dist_model()
# register_nearest_neighbor_dist_model()

cat("Model info for 'nearest_neighbor_dist':\n")
parsnip::show_model_info("nearest_neighbor_dist")
cat("\n")

## -----------------------------------------------------------------------
## 3. Simple classification dataset (mixed or numeric, here: iris)
## -----------------------------------------------------------------------

dat <- as_tibble(iris) |>
  mutate(Species = factor(Species))

set.seed(123)
split <- initial_split(dat, prop = 0.75, strata = Species)
train <- training(split)
test  <- testing(split)

cat("Train n =", nrow(train), " | Test n =", nrow(test), "\n\n")

## -----------------------------------------------------------------------
## 4. Recipe: standard tidymodels formula
## -----------------------------------------------------------------------

rec <- recipe(Species ~ ., data = train)

## -----------------------------------------------------------------------
## 5. Model spec: nearest_neighbor_dist + mdist as distance
## -----------------------------------------------------------------------
## IMPORTANT:
##  - manydist::mdist must accept (train, new_data = , ...)
##  - and return an MDist object with $distance (n_test x n_train)

mdl <- nearest_neighbor_dist(
  mode      = "classification",
  neighbors = 5,
  dist_fun  = mdist,         # dev mdist from your package
  dist_args = list(preset = "gower")   # or your preferred preset
) |>
  set_engine("precomputed")

print(mdl)

## -----------------------------------------------------------------------
## 6. Workflow: recipe + model
## -----------------------------------------------------------------------

wf <- workflow() |>
  add_recipe(rec) |>
  add_model(mdl)

print(wf)

## -----------------------------------------------------------------------
## 7. Fit the workflow
## -----------------------------------------------------------------------

fit_wf <- fit(wf, data = train)

cat("\nFit completed. Workflow fit object:\n")
print(fit_wf)

## -----------------------------------------------------------------------
## 8. Predictions and simple accuracy check
## -----------------------------------------------------------------------

# class predictions
pred_class <- predict(fit_wf, new_data = test, type = "class") |>
  bind_cols(test |> select(Species))

# prob predictions
pred_prob <- predict(fit_wf, new_data = test, type = "prob")

cat("\nHead of class predictions:\n")
print(head(pred_class))

cat("\nHead of probability predictions:\n")
print(head(pred_prob))

acc <- yardstick::accuracy(pred_class, truth = Species, estimate = .pred_class)
cat("\nClassification accuracy (nearest_neighbor_dist + mdist):\n")
print(acc)
cat("\n")

## -----------------------------------------------------------------------
## 9. Optional: quick sanity check with yardstick confusion matrix
## -----------------------------------------------------------------------

cat("Confusion matrix:\n")
print(yardstick::conf_mat(pred_class, truth = Species, estimate = .pred_class))
cat("\nAll good if you got accuracy and a confusion matrix without errors.\n")
