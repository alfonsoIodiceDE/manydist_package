## =====================================================================
##  CLEAN TEST: nearest_neighbor_dist + mdist + tidymodels
##  (All parsnip glue defined *here*, no package registration used)
## =====================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidymodels)
  library(pkgload)
  library(dplyr)
  library(tibble)
})

## -----------------------------------------------------------------------
## 1. Load dev manydist (ONLY)
## -----------------------------------------------------------------------

pkg_path <- "../manydist"

if ("manydist" %in% loadedNamespaces()) {
  unloadNamespace("manydist")
}

pkgload::load_all(pkg_path, export_all = TRUE, helpers = TRUE, attach_testthat = FALSE)

cat("Loaded manydist from:\n  ",
    getNamespaceInfo(asNamespace("manydist"), "path"), "\n\n", sep = "")

# We’ll use mdist from the package, but NOT the package’s registration


## -----------------------------------------------------------------------
## 2. Define nearest_neighbor_dist model spec + helpers LOCALLY
##    (These live in .GlobalEnv for this test)
## -----------------------------------------------------------------------

nearest_neighbor_dist <- function(mode = "classification",
                                  neighbors = NULL,
                                  dist_fun  = NULL,
                                  dist_args = NULL) {
  parsnip::new_model_spec(
    cls  = "nearest_neighbor_dist",
    args = list(
      neighbors = rlang::enquo(neighbors),
      dist_fun  = rlang::enquo(dist_fun),
      dist_args = rlang::enquo(dist_args)
    ),
    eng_args = NULL,
    mode  = mode,
    method = NULL,
    engine = NULL
  )
}

# --- kNN helpers on distance matrix ------------------------------------

.knn_class_from_dist <- function(D, y_train, k) {
  nn_idx <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  lv <- levels(y_train)
  out <- apply(nn_idx, 2L, function(idx) {
    tab <- table(y_train[idx])
    names(tab)[which.max(tab)]
  })
  factor(out, levels = lv)
}

.knn_prob_from_dist <- function(D, y_train, k) {
  classes <- levels(y_train)
  nn_idx  <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  probs <- t(apply(nn_idx, 2L, function(idx) {
    tab <- table(factor(y_train[idx], levels = classes))
    as.numeric(tab / sum(tab))
  }))
  colnames(probs) <- classes
  as.data.frame(probs)
}

.knn_reg_from_dist <- function(D, y_train, k) {
  nn_idx <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  apply(nn_idx, 2L, function(idx) mean(y_train[idx]))
}

.dist_to_matrix <- function(D_raw) {
  if (inherits(D_raw, "MDist")) {
    as.matrix(D_raw$distance)
  } else if (inherits(D_raw, "dist")) {
    as.matrix(D_raw)
  } else if (is.matrix(D_raw)) {
    D_raw
  } else {
    stop("dist_fun must return an 'MDist', 'dist', or numeric matrix.",
         call. = FALSE)
  }
}

# --- fit + predict (parsnip-facing) ------------------------------------

fit_nearest_neighbor_dist <- function(x, y,
                                      neighbors = 5,
                                      dist_fun,
                                      dist_args = list()) {

  if (rlang::is_quosure(neighbors)) neighbors <- rlang::eval_tidy(neighbors)
  if (rlang::is_quosure(dist_fun))  dist_fun  <- rlang::eval_tidy(dist_fun)
  if (rlang::is_quosure(dist_args)) dist_args <- rlang::eval_tidy(dist_args)

  if (missing(dist_fun) || is.null(dist_fun)) {
    stop("Please supply 'dist_fun(train, new_data, ...)' returning an MDist or matrix.",
         call. = FALSE)
  }

  y_cols <- character(0)
  if (is.data.frame(y)) {
    y_cols <- names(y)
    y      <- y[[1L]]
  }

  if (is.data.frame(x) && length(y_cols)) {
    x <- x[, setdiff(names(x), y_cols), drop = FALSE]
  }

  obj <- list(
    x         = x,
    y         = y,
    k         = neighbors,
    dist_fun  = dist_fun,
    dist_args = dist_args
  )
  attr(obj, "y_cols") <- y_cols

  structure(obj, class = "nearest_neighbor_dist")
}

predict_nearest_neighbor_dist_class <- function(object,
                                                new_data,
                                                type = c("class", "prob")) {
  type <- match.arg(type)

  y_cols <- attr(object, "y_cols", exact = TRUE)
  if (is.data.frame(new_data) && !is.null(y_cols)) {
    new_data <- new_data[, setdiff(names(new_data), y_cols), drop = FALSE]
  }

  D_raw <- do.call(object$dist_fun,
                   c(list(object$x, new_data = new_data), object$dist_args))
  D <- .dist_to_matrix(D_raw)

  if (type == "class") {
    .knn_class_from_dist(D, object$y, object$k)
  } else {
    .knn_prob_from_dist(D, object$y, object$k)
  }
}

predict_nearest_neighbor_dist_reg <- function(object, new_data) {
  y_cols <- attr(object, "y_cols", exact = TRUE)
  if (is.data.frame(new_data) && !is.null(y_cols)) {
    new_data <- new_data[, setdiff(names(new_data), y_cols), drop = FALSE]
  }

  D_raw <- do.call(object$dist_fun,
                   c(list(object$x, new_data = new_data), object$dist_args))
  D <- .dist_to_matrix(D_raw)

  .knn_reg_from_dist(D, object$y, object$k)
}

## -----------------------------------------------------------------------
## 3. DEV registration (pkg = NULL → global env)
## -----------------------------------------------------------------------

register_nearest_neighbor_dist_model_dev <- function() {
  q <- function(expr) try(suppressWarnings(expr), silent = TRUE)

  # model + modes
  q(parsnip::set_new_model("nearest_neighbor_dist"))
  q(parsnip::set_model_mode("nearest_neighbor_dist", "classification"))
  q(parsnip::set_model_mode("nearest_neighbor_dist", "regression"))

  # engine
  q(parsnip::set_model_engine("nearest_neighbor_dist", "classification", "precomputed"))
  q(parsnip::set_model_engine("nearest_neighbor_dist", "regression",     "precomputed"))

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

  # fit (pkg = NULL → look in .GlobalEnv)
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

  # predict
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

register_nearest_neighbor_dist_model_dev()
parsnip::show_model_info("nearest_neighbor_dist")

## -----------------------------------------------------------------------
## 4. Test on iris (classification)
## -----------------------------------------------------------------------

dat <- as_tibble(iris) |>
  mutate(Species = factor(Species))

set.seed(123)
split <- initial_split(dat, prop = 0.75, strata = Species)
train <- training(split)
test  <- testing(split)

rec <- recipe(Species ~ ., data = train)

mdl <- nearest_neighbor_dist(
  mode      = "classification",
  neighbors = 5,
  dist_fun  = mdist,                # from dev manydist
  dist_args = list(preset = "gower")
) |>
  set_engine("precomputed")

wf <- workflow() |>
  add_recipe(rec) |>
  add_model(mdl)

fit_wf <- fit(wf, data = train)

pred_class <- predict(fit_wf, new_data = test, type = "class") |>
  bind_cols(test |> select(Species))

pred_prob <- predict(fit_wf, new_data = test, type = "prob")

cat("\nHead of class predictions:\n")
print(head(pred_class))

cat("\nHead of probability predictions:\n")
print(head(pred_prob))

acc <- yardstick::accuracy(pred_class, truth = Species, estimate = .pred_class)
cat("\nAccuracy:\n")
print(acc)

cat("\nConfusion matrix:\n")
print(yardstick::conf_mat(pred_class, truth = Species, estimate = .pred_class))
