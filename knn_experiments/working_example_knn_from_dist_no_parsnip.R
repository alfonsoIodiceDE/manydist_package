## --- Core kNN helpers on a distance matrix ---------------------------------

# majority vote (classification)
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

# posterior probabilities
.knn_prob_from_dist <- function(D, y_train, k) {
  classes <- levels(y_train)
  nn_idx <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  probs <- t(apply(nn_idx, 2L, function(idx) {
    tab <- table(factor(y_train[idx], levels = classes))
    as.numeric(tab / sum(tab))
  }))
  colnames(probs) <- classes
  as.data.frame(probs)
}

# regression mean
.knn_reg_from_dist <- function(D, y_train, k) {
  nn_idx <- apply(D, 1L, function(d) order(d)[seq_len(k)])
  if (is.vector(nn_idx)) nn_idx <- matrix(nn_idx, nrow = 1L)
  apply(nn_idx, 2L, function(idx) mean(y_train[idx]))
}

# minimal “fit” object: we just store train + distance machinery
fit_knn_dist <- function(x, y, k = 5, dist_fun, dist_args = list()) {
  stopifnot(is.function(dist_fun), k >= 1)
  structure(
    list(
      x         = x,
      y         = y,
      k         = k,
      dist_fun  = dist_fun,
      dist_args = dist_args
    ),
    class = "knn_dist"
  )
}

predict_knn_dist_class <- function(object, new_data, type = c("class","prob")) {
  type <- match.arg(type)
  D_raw <- do.call(object$dist_fun,
                   c(list(object$x, new_data = new_data), object$dist_args))

  # Coerce MDist → matrix
  if (inherits(D_raw, "MDist")) {
    D <- as.matrix(D_raw$distance)
  } else if (is.matrix(D_raw)) {
    D <- D_raw
  } else {
    stop("dist_fun must return an MDist object or a numeric matrix.")
  }

  if (type == "class") .knn_class_from_dist(D, object$y, object$k)
  else                 .knn_prob_from_dist(D, object$y, object$k)
}

predict_knn_dist_reg <- function(object, new_data) {
  D_raw <- do.call(object$dist_fun,
                   c(list(object$x, new_data = new_data), object$dist_args))

  if (inherits(D_raw, "MDist")) {
    D <- as.matrix(D_raw$distance)
  } else if (is.matrix(D_raw)) {
    D <- D_raw
  } else {
    stop("dist_fun must return an MDist object or a numeric matrix.")
  }

  .knn_reg_from_dist(D, object$y, object$k)
}

library(dplyr)
library(rsample)
library(yardstick)
library(tibble)

dummy_dist <- function(train, new_data, ...) {
  X <- as.matrix(dplyr::select(train, -Species))
  Z <- as.matrix(dplyr::select(new_data, -Species))
  Dfull <- as.matrix(dist(rbind(Z, X)))  # (n_test + n_train) × (n_test + n_train)
  Dfull[seq_len(nrow(Z)), (nrow(Z) + 1):(nrow(Z) + nrow(X))]  # n_test × n_train
}


dat <- as_tibble(iris) %>% mutate(Species = factor(Species))
set.seed(1)
sp  <- initial_split(dat, prop = 0.7, strata = Species)
tr  <- training(sp)
te  <- testing(sp)

fit0 <- fit_knn_dist(x = tr, y = tr$Species, k = 5,
                     dist_fun = dummy_dist, dist_args = list())

pred_class <- predict_knn_dist_class(fit0, te, type = "class")
acc <- accuracy_vec(truth = te$Species, estimate = pred_class)
cat("Dummy kNN accuracy:", acc, "\n")

# try to do the same with mdist
source("R/mdist.R")

fit_m <- fit_knn_dist(
  x         = tr,
  y         = tr$Species,
  k         = 5,
  dist_fun  = mdist,
  dist_args = list(preset = "euclidean_onehot")  # or whatever your mdist takes
)

pred_class_m <- predict_knn_dist_class(fit_m, te, type = "prob")


nearest_neighbor_dist <- function(mode      = "classification",
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

  parsnip::new_model_spec(
    "nearest_neighbor_dist",
    args     = args,
    eng_args = NULL,
    mode     = mode,
    method   = NULL,
    engine   = NULL
  )
}


dummy_dist <- function(train_x, new_data, ...) {
  X <- as.matrix(train_x)
  Z <- as.matrix(new_data)
  Dfull <- as.matrix(dist(rbind(Z, X)))  # (n_test + n_train) × (n_test + n_train)
  Dfull[seq_len(nrow(Z)), (nrow(Z) + 1):(nrow(Z) + nrow(X))]  # n_test × n_train
}


register_nearest_neighbor_dist_dev <- function() {
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

  # arguments
  q(parsnip::set_model_arg(
    model    = "nearest_neighbor_dist",
    eng      = "precomputed",
    parsnip  = "neighbors",
    original = "k",   # <- argument name in fit_knn_dist()
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

  # FIT: use your existing fit_knn_dist()
  q(parsnip::set_fit(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "classification",
    value = list(
      interface = "data.frame",
      protect   = c("x", "y"),
      func      = c(pkg = NULL, fun = "fit_knn_dist"),
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
      func      = c(pkg = NULL, fun = "fit_knn_dist"),
      defaults  = list()
    )
  ))

  # PREDICT: classification
  q(parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "classification",
    type  = "class",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = NULL, fun = "predict_knn_dist_class"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
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
      pre  = NULL,
      post = NULL,
      func = c(pkg = NULL, fun = "predict_knn_dist_class"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
        type     = "prob"
      )
    )
  ))

  # PREDICT: regression
  q(parsnip::set_pred(
    model = "nearest_neighbor_dist",
    eng   = "precomputed",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = NULL, fun = "predict_knn_dist_reg"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  ))

  invisible(TRUE)
}

# call it
register_nearest_neighbor_dist_dev()
parsnip::show_model_info("nearest_neighbor_dist")


library(tidymodels)

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
  dist_args = list()
) |>
  set_engine("precomputed")

# Formula interface: parsnip will pass x = predictors, y = Species
fit_parsnip <- fit(
  spec,
  Species ~ .,
  data = tr
)

