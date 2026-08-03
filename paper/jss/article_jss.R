options(prompt = "R> ", continue = "+  ", width = 70,
        useFancyQuotes = FALSE)


###################################################
### Section 2: A unified mixed-variable distance framework
###################################################

penguins <- palmerpenguins::penguins |>
  dplyr::select(
    species,
    bill_length_mm,
    bill_depth_mm,
    flipper_length_mm,
    body_mass_g,
    island,
    sex
  ) |>
  tidyr::drop_na()

penguins_x <- penguins |>
  dplyr::select(-species)
penguins_y <- penguins |>
  dplyr::pull(species)

dependence_data <- penguins |>
  dplyr::transmute(
    penguin_id = dplyr::row_number(),
    flipper_length_mm,
    body_mass_g,
    flipper_z = as.numeric(scale(.data$flipper_length_mm)),
    body_mass_z = as.numeric(scale(.data$body_mass_g))
  )

standardized_data <- dependence_data |>
  dplyr::select(.data$flipper_z, .data$body_mass_z) |>
  as.matrix()

pca_fit <- stats::prcomp(
  standardized_data,
  center = FALSE,
  scale. = FALSE
)

whitened_scores <- sweep(
  pca_fit$x[, 1:2, drop = FALSE],
  MARGIN = 2,
  STATS = pca_fit$sdev[1:2],
  FUN = "/"
)

dependence_data <- dependence_data |>
  dplyr::mutate(
    pc1 = pca_fit$x[, 1],
    pc2 = pca_fit$x[, 2],
    whitened_pc1 = whitened_scores[, 1],
    whitened_pc2 = whitened_scores[, 2]
  )

pc1_targets <- stats::quantile(
  dependence_data$pc1,
  probs = c(0.2, 0.8)
)

selected_ids <- purrr::map_int(pc1_targets, function(target) {
  dependence_data |>
    dplyr::mutate(target_distance = abs(.data$pc1 - target)) |>
    dplyr::slice_min(
      order_by = .data$target_distance,
      n = 20,
      with_ties = FALSE
    ) |>
    dplyr::slice_min(
      order_by = abs(.data$pc2),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::pull(.data$penguin_id)
})

selected_points <- dependence_data |>
  dplyr::filter(.data$penguin_id %in% selected_ids) |>
  dplyr::arrange(.data$flipper_z) |>
  dplyr::mutate(point = c("A", "B"))

raw_segment <- selected_points |>
  dplyr::summarise(
    x = dplyr::first(.data$flipper_length_mm),
    y = dplyr::first(.data$body_mass_g),
    xend = dplyr::last(.data$flipper_length_mm),
    yend = dplyr::last(.data$body_mass_g),
    distance = sqrt((.data$xend - .data$x)^2 + (.data$yend - .data$y)^2),
    xmid = (.data$x + .data$xend) / 2,
    ymid = (.data$y + .data$yend) / 2
  )

standardized_segment <- selected_points |>
  dplyr::summarise(
    x = dplyr::first(.data$flipper_z),
    y = dplyr::first(.data$body_mass_z),
    xend = dplyr::last(.data$flipper_z),
    yend = dplyr::last(.data$body_mass_z),
    distance = sqrt((.data$xend - .data$x)^2 + (.data$yend - .data$y)^2),
    xmid = (.data$x + .data$xend) / 2,
    ymid = (.data$y + .data$yend) / 2
  )

component_segments <- tibble::tibble(
  x = c(standardized_segment$x, standardized_segment$xend),
  y = c(standardized_segment$y, standardized_segment$y),
  xend = c(standardized_segment$xend, standardized_segment$xend),
  yend = c(standardized_segment$y, standardized_segment$yend)
)

mahalanobis_segment <- selected_points |>
  dplyr::summarise(
    x = dplyr::first(.data$whitened_pc1),
    y = dplyr::first(.data$whitened_pc2),
    xend = dplyr::last(.data$whitened_pc1),
    yend = dplyr::last(.data$whitened_pc2),
    distance = sqrt((.data$xend - .data$x)^2 + (.data$yend - .data$y)^2),
    xmid = (.data$x + .data$xend) / 2,
    ymid = (.data$y + .data$yend) / 2
  )

common_score_limits <- c(-3.5, 3.5)
common_score_breaks <- seq(-3, 3, by = 1)

distance_plot_theme <- ggplot2::theme_minimal(base_size = 8) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 8),
    axis.text = ggplot2::element_text(size = 7),
    plot.margin = ggplot2::margin(3, 4, 3, 3)
  )

raw_plot <- ggplot2::ggplot(
  dependence_data,
  ggplot2::aes(x = .data$flipper_length_mm, y = .data$body_mass_g)
) +
  ggplot2::geom_point(size = 1.1, color = "grey35", alpha = 0.55) +
  ggplot2::stat_ellipse(
    type = "norm",
    level = 0.68,
    linewidth = 0.5,
    color = "#A33F3F"
  ) +
  ggplot2::geom_segment(
    data = raw_segment,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    inherit.aes = FALSE,
    linewidth = 0.65,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"))
  ) +
  ggplot2::geom_point(
    data = selected_points,
    size = 2.1,
    color = "#A33F3F"
  ) +
  ggplot2::geom_text(
    data = selected_points,
    ggplot2::aes(label = .data$point),
    nudge_y = 170,
    size = 2.7,
    fontface = "bold"
  ) +
  ggplot2::geom_label(
    data = raw_segment,
    ggplot2::aes(
      x = .data$xmid,
      y = .data$ymid,
      label = paste0("d[E]^{raw} == ", round(.data$distance, 0))
    ),
    parse = TRUE,
    inherit.aes = FALSE,
    size = 2.3,
    linewidth = 0,
    fill = "white"
  ) +
  ggplot2::labs(
    x = "Flipper length (mm)",
    y = "Body mass (g)"
  ) +
  distance_plot_theme

independence_plot <- ggplot2::ggplot(
  dependence_data,
  ggplot2::aes(x = .data$flipper_z, y = .data$body_mass_z)
) +
  ggplot2::geom_point(size = 0.65, color = "grey70") +
  ggplot2::stat_ellipse(
    type = "norm",
    level = 0.68,
    linewidth = 0.5,
    color = "#A33F3F"
  ) +
  ggplot2::geom_segment(
    data = component_segments,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    inherit.aes = FALSE,
    linewidth = 0.55,
    linetype = "dashed",
    color = "#A33F3F"
  ) +
  ggplot2::geom_segment(
    data = standardized_segment,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    inherit.aes = FALSE,
    linewidth = 0.65,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"))
  ) +
  ggplot2::geom_point(
    data = selected_points,
    size = 1.7,
    color = "#A33F3F"
  ) +
  ggplot2::geom_text(
    data = selected_points,
    ggplot2::aes(label = .data$point),
    nudge_y = 0.2,
    size = 2.7,
    fontface = "bold"
  ) +
  ggplot2::geom_label(
    data = standardized_segment,
    ggplot2::aes(
      x = .data$xmid,
      y = .data$ymid,
      label = paste0("d[E]^{z} == ", round(.data$distance, 2))
    ),
    parse = TRUE,
    inherit.aes = FALSE,
    size = 2.3,
    linewidth = 0,
    fill = "white"
  ) +
  ggplot2::scale_x_continuous(
    limits = common_score_limits,
    breaks = common_score_breaks,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::scale_y_continuous(
    limits = common_score_limits,
    breaks = common_score_breaks,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Flipper length (z-score)",
    y = "Body mass (z-score)"
  ) +
  distance_plot_theme

association_plot <- ggplot2::ggplot(
  dependence_data,
  ggplot2::aes(x = .data$whitened_pc1, y = .data$whitened_pc2)
) +
  ggplot2::geom_point(size = 0.65, color = "grey70") +
  ggplot2::stat_ellipse(
    type = "norm",
    level = 0.68,
    linewidth = 0.5,
    color = "#2A7F62"
  ) +
  ggplot2::geom_segment(
    data = mahalanobis_segment,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    inherit.aes = FALSE,
    linewidth = 0.65,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"))
  ) +
  ggplot2::geom_point(
    data = selected_points,
    size = 1.7,
    color = "#2A7F62"
  ) +
  ggplot2::geom_text(
    data = selected_points,
    ggplot2::aes(label = .data$point),
    nudge_y = 0.2,
    size = 2.7,
    fontface = "bold"
  ) +
  ggplot2::geom_label(
    data = mahalanobis_segment,
    ggplot2::aes(
      x = .data$xmid,
      y = .data$ymid,
      label = paste0("d[M] == ", round(.data$distance, 2))
    ),
    parse = TRUE,
    inherit.aes = FALSE,
    size = 2.3,
    linewidth = 0,
    fill = "white"
  ) +
  ggplot2::scale_x_continuous(
    limits = common_score_limits,
    breaks = common_score_breaks,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::scale_y_continuous(
    limits = common_score_limits,
    breaks = common_score_breaks,
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Whitened PC 1",
    y = "Whitened PC 2"
  ) +
  distance_plot_theme

raw_plot
independence_plot
association_plot

###################################################
### Section 3: Distance construction
###################################################

library("manydist")

d_u_dep <- mdist(
  penguins_x,
  preset = "u_dep"
)

d_u_dep

as.matrix(d_u_dep$to_dist())[1:4, 1:4] |>
  round(digits = 2)

d_custom <- mdist(
  penguins_x,
  preset = "custom",
  method_cat = "matching",
  method_num = "pc_scores",
  commensurable = TRUE
)

## Commensurability cancels any per-variable affine scaling, so the sd,
## range and robust settings of method_num return the same distance when
## commensurable = TRUE. Both comparisons below should return TRUE.

purrr::map(
  c("sd", "range", "robust"),
  ~ mdist(
    penguins_x,
    preset = "custom",
    method_num = .x,
    method_cat = "matching",
    commensurable = TRUE
  )$to_dist()
) |>
  (\(d) c(
    all.equal(d[[1]], d[[2]]),
    all.equal(d[[1]], d[[3]])
  ))()

## The same three settings do change the distance when the weights are
## not applied. Contrast with the check above.

purrr::map(
  c("sd", "range", "robust"),
  ~ mdist(
    penguins_x,
    preset = "custom",
    method_num = .x,
    method_cat = "matching",
    commensurable = FALSE
  )$to_dist()
) |>
  (\(d) c(
    isTRUE(all.equal(d[[1]], d[[2]])),
    isTRUE(all.equal(d[[1]], d[[3]]))
  ))()

d_response <- mdist(
  penguins,
  response = species,
  preset = "u_dep"
)

###################################################
### Section 4: Distance diagnostics and benchmarking
###################################################

lovo_gower <- lovo_mdist(
  penguins,
  response = species,
  response_used = FALSE,
  preset = "gower"
)

lovo_gower$results |>
  dplyr::select(
    variable,
    relative_distance,
    ac_importance
  ) |>
  dplyr::arrange(dplyr::desc(relative_distance))

lovo_gower$autoplot(
  metric = "relative_distance",
  reorder = TRUE
) +
  ggplot2::guides(colour = "none")

lovo_comparison <- compare_lovo_mdist(
  x = penguins_x,
  methods = list(
    Gower = list(preset = "gower"),
    Euclidean = list(preset = "euclidean"),
    `Unbiased independent` = list(preset = "u_indep"),
    `Unbiased dependent` = list(preset = "u_dep")
  ),
  cluster_k = 3,
  cluster_methods = "pam"
)

lovo_comparison$autoplot(
  metric = "ari_pam",
  reorder = TRUE
) +
  ggplot2::labs(
    x = NULL,
    title = NULL,
    color = "Distance"
  ) +
  ggplot2::theme(
    legend.position = "top"
  )

candidate_specs <- all_dist_method_specs(
  mode = "presets_only",
  preset = c("gower", "u_indep", "u_dep")
) |>
  dplyr::mutate(
    label = dplyr::recode_values(
      preset,
      "gower" ~ "Gower",
      "u_indep" ~ "Unbiased independent",
      "u_dep" ~ "Unbiased dependent",
      unmatched = "error"
    )
  )

set.seed(2026)

distance_benchmark <- benchmark_mdist(
  x = penguins_x,
  specs = candidate_specs,
  cluster_k = 3,
  cluster_methods = "pam"
)

benchmark_pairs <- benchmark_comparisons(distance_benchmark) |>
  dplyr::select(
    method_1,
    method_2,
    mad,
    relative_distance,
    alienation,
    ari_pam
  )

benchmark_pairs |>
  dplyr::transmute(
    Pair = paste(method_1, method_2, sep = " -- "),
    MAD = mad,
    `Relative distance` = relative_distance,
    Alienation = alienation,
    `PAM ARI` = ari_pam
  ) |>
  knitr::kable(
    digits = 3,
    align = c("l", "r", "r", "r", "r")
  )

## Why the two Gower rows of the table above share a MAD and a relative
## distance. Commensurability sets the expected contribution of every
## variable to a common constant, so the mean pairwise distance of any
## commensurable specification on these six predictors is the same. The
## two u_ entries below should therefore be identical, and their common
## value minus the Gower mean reproduces the shared MAD. Raw magnitude
## against a non-commensurable baseline cannot separate two
## commensurable distances; only the geometry measures can.

mean_pairwise <- c(
  u_indep = mean(mdist(penguins_x, preset = "u_indep")$to_dist()),
  u_dep = mean(mdist(penguins_x, preset = "u_dep")$to_dist()),
  gower = mean(mdist(penguins_x, preset = "gower")$to_dist())
)

round(mean_pairwise, 3)

round(
  c(
    u_indep_minus_gower = unname(mean_pairwise["u_indep"] -
                                   mean_pairwise["gower"]),
    u_dep_minus_gower = unname(mean_pairwise["u_dep"] -
                                 mean_pairwise["gower"])
  ),
  3
)

###################################################
### Section 5: Distance-based learning pipelines
###################################################

data("wdi_2022", package = "manydist")

wdi_clustering <- wdi_2022 |>
  dplyr::select(-country)

distance_presets <- c(
  gower = "gower",
  euclidean = "euclidean",
  u_indep = "u_indep",
  u_dep = "u_dep",
  u_mix = "u_mix"
)

distance_labels <- c(
  gower = "Gower",
  euclidean = "Euclidean",
  u_indep = "Unbiased independent",
  u_dep = "Unbiased dependent",
  u_mix = "Unbiased mixed"
)

## The neighbourhood grid is defined here because both the refitting
## comparison below and the supervised tuning further down use it.

neighbor_grid <- tibble::tibble(
  neighbors = c(1L, 3L, 5L, 7L, 9L, 15L)
)

### Refitting the distance within resamples

## The classification data and the resampling scheme are prepared first
## because both the refitting comparison below and the supervised task
## further down use them.

wdi_classification <- wdi_2022 |>
  dplyr::filter(income != "Not classified") |>
  dplyr::mutate(
    income = forcats::fct_drop(income),
    region = forcats::fct_collapse(
      region,
      Americas = c(
        "North America",
        "Latin America & Caribbean"
      )
    )
  ) |>
  dplyr::select(-country, -lending)

wdi_classification |>
  dplyr::count(income)

set.seed(2026)

wdi_split <- rsample::initial_split(
  wdi_classification,
  prop = 0.75,
  strata = income
)

wdi_train <- rsample::training(wdi_split)

set.seed(2026)

wdi_folds <- rsample::vfold_cv(
  wdi_train,
  v = 5,
  strata = income
)

## Three cross-validated estimates of the same nearest-neighbour
## classifier. leaky_cv() and honest_cv() share the voting code in
## knn_vote() and differ ONLY in whether the distance is estimated once
## on the whole training partition or refit on each analysis set, so
## their difference isolates the leakage. refit_cv() runs the
## step_mdist() pipeline and checks that honest_cv() reproduces what the
## package does internally.
##
## fit_mdist() drops `income` from the predictors in every case and
## supplies it as a tidyselect `response` only for the presets that use
## it. Passing the full frame without dropping it would silently make
## the outcome an ordinary matching predictor.

response_aware_presets <- c("u_dep", "u_mix")
leakage_presets <- c("u_indep", "u_mix", "u_dep")

knn_vote <- function(d_block, train_y, k) {
  apply(d_block, 1, function(row) {
    nn <- order(row)[seq_len(k)]
    names(sort(table(train_y[nn]), decreasing = TRUE))[1]
  })
}

fit_mdist <- function(data, preset, new_data = NULL) {
  args <- if (preset %in% response_aware_presets) {
    list(data, response = quote(income))
  } else {
    list(dplyr::select(data, -income))
  }
  args$preset <- preset
  if (!is.null(new_data)) {
    args$new_data <- new_data
  }
  do.call(mdist, args)
}

## Assertion. All three presets must see the same predictor set, with
## `income` absent from every one of them.

purrr::map(leakage_presets, ~ names(fit_mdist(wdi_train, .x)$data))

## Assertion. With `region` the only categorical predictor, the
## association-based presets can only condition on the outcome.
## u_indep should report "matching"; u_dep and u_mix should report
## "tvd".

purrr::map_chr(
  leakage_presets,
  ~ fit_mdist(wdi_train, .x)$params$method_cat
)

fold_accuracy <- function(d_block, train_y, truth, k) {
  pred <- knn_vote(d_block, train_y, k)
  yardstick::bal_accuracy_vec(
    truth = truth,
    estimate = factor(pred, levels = levels(wdi_train$income))
  )
}

leaky_cv <- function(preset, k) {
  m_all <- as.matrix(fit_mdist(wdi_train, preset)$to_dist())

  purrr::map_dbl(wdi_folds$splits, function(split) {
    fit_id <- split$in_id
    out_id <- rsample::complement(split)
    fold_accuracy(
      m_all[out_id, fit_id, drop = FALSE],
      wdi_train$income[fit_id],
      wdi_train$income[out_id],
      k
    )
  })
}

honest_cv <- function(preset, k) {
  purrr::map_dbl(wdi_folds$splits, function(split) {
    fit_data <- rsample::analysis(split)
    out_data <- rsample::assessment(split)

    d_block <- unclass(
      fit_mdist(
        fit_data,
        preset,
        new_data = dplyr::select(out_data, -income)
      )$distance
    )

    stopifnot(
      nrow(d_block) == nrow(out_data),
      ncol(d_block) == nrow(fit_data)
    )

    fold_accuracy(
      d_block,
      fit_data$income,
      out_data$income,
      k
    )
  })
}

refit_cv <- function(preset, k) {
  wf <- workflows::workflow() |>
    workflows::add_recipe(
      recipes::recipe(income ~ ., data = wdi_train) |>
        step_mdist(
          recipes::all_predictors(),
          preset = preset,
          output = "distance_to_training"
        )
    ) |>
    workflows::add_model(
      nearest_neighbor_dist(
        mode = "classification",
        neighbors = k
      ) |>
        parsnip::set_engine("manydist")
    )

  tune::fit_resamples(
    wf,
    resamples = wdi_folds,
    metrics = yardstick::metric_set(yardstick::bal_accuracy)
  ) |>
    tune::collect_metrics() |>
    dplyr::filter(.metric == "bal_accuracy") |>
    dplyr::pull(mean)
}

set.seed(2026)

leakage_folds <- tibble::tibble(preset = leakage_presets) |>
  dplyr::mutate(
    leaky = purrr::map(preset, leaky_cv, k = 5L),
    honest = purrr::map(preset, honest_cv, k = 5L),
    pipeline = purrr::map_dbl(preset, refit_cv, k = 5L)
  )

## The folds are paired, so the standard error of the fold-wise
## differences is the relevant one, not the standard error of either
## estimate on its own.

leakage_comparison <- leakage_folds |>
  dplyr::mutate(
    single_matrix = purrr::map_dbl(leaky, mean),
    refit_per_fold = purrr::map_dbl(honest, mean),
    optimism = single_matrix - refit_per_fold,
    paired_se = purrr::map2_dbl(
      leaky,
      honest,
      ~ stats::sd(.x - .y) / sqrt(length(.x))
    ),
    pipeline_gap = refit_per_fold - pipeline
  )

knitr::kable(
  leakage_comparison |>
    dplyr::transmute(
      Distance = unname(distance_labels[preset]),
      `Single matrix` = single_matrix,
      `Refit per fold` = refit_per_fold,
      Optimism = optimism,
      `Paired SE` = paired_se
    ),
  digits = 3
)

## Check. honest_cv() should reproduce the step_mdist() pipeline, so
## pipeline_gap should be at or near zero. A non-zero value is the
## difference between knn_vote() and the packaged engine, and would have
## to be subtracted from the optimism column above.

leakage_comparison |>
  dplyr::transmute(
    Distance = unname(distance_labels[preset]),
    honest = refit_per_fold,
    pipeline,
    gap = pipeline_gap
  )

## Check. Fold-by-fold differences, to see whether the optimism is
## consistent in sign across resamples or driven by a single fold.

leakage_folds |>
  dplyr::mutate(
    fold = list(seq_len(nrow(wdi_folds))),
    difference = purrr::map2(leaky, honest, ~ .x - .y)
  ) |>
  dplyr::select(preset, fold, difference) |>
  tidyr::unnest(c(fold, difference)) |>
  tidyr::pivot_wider(
    names_from = preset,
    values_from = difference
  )

## Does the leak change which specification would be selected? Repeat
## the comparison across the whole neighbourhood grid and record the
## rank each specification receives under either approach.

set.seed(2026)

leakage_grid <- tidyr::expand_grid(
  preset = leakage_presets,
  k = neighbor_grid$neighbors
) |>
  dplyr::mutate(
    single_matrix = purrr::map2_dbl(preset, k, ~ mean(leaky_cv(.x, .y))),
    refit_per_fold = purrr::map2_dbl(preset, k, ~ mean(honest_cv(.x, .y)))
  ) |>
  dplyr::group_by(k) |>
  dplyr::mutate(
    leaky_rank = rank(-single_matrix),
    honest_rank = rank(-refit_per_fold)
  ) |>
  dplyr::ungroup()

leakage_grid |>
  dplyr::transmute(
    Distance = unname(distance_labels[preset]),
    Neighbours = k,
    `Single matrix` = leaky_rank,
    `Refit per fold` = honest_rank
  ) |>
  tidyr::pivot_wider(
    names_from = Neighbours,
    values_from = c(`Single matrix`, `Refit per fold`),
    names_sep = ", k = "
  ) |>
  knitr::kable()

## The neighbourhood sizes at which the two approaches would select a
## different specification.

leakage_grid |>
  dplyr::group_by(k) |>
  dplyr::summarise(
    leaky_best = unname(distance_labels[preset[which.max(single_matrix)]]),
    honest_best = unname(distance_labels[preset[which.max(refit_per_fold)]]),
    .groups = "drop"
  ) |>
  dplyr::mutate(disagree = leaky_best != honest_best)

### Unsupervised country clustering

clustering_recipes <- purrr::map(
  distance_presets,
  ~ recipes::recipe(~ ., data = wdi_clustering) |>
    step_mdist(
      recipes::all_predictors(),
      preset = .x,
      output = "pairwise"
    )
)

wdi_pam_fit <- generics::fit(
  pam_dist(num_clusters = 4),
  recipe = clustering_recipes$u_mix,
  data = wdi_clustering
)

wdi_spectral_fit <- generics::fit(
  spectral_dist(num_clusters = 4, nstart = 50),
  recipe = clustering_recipes$u_mix,
  data = wdi_clustering
)

fit_clustering_pipeline <- function(distance, method) {
  model <- if (method == "PAM") {
    pam_dist(num_clusters = 4)
  } else {
    spectral_dist(num_clusters = 4, nstart = 50)
  }

  fitted <- generics::fit(
    model,
    recipe = clustering_recipes[[distance]],
    data = wdi_clustering
  )

  assignment <- predict(fitted)$.pred_cluster |>
    as.integer()

  silhouette <- cluster::silhouette(
    assignment,
    stats::as.dist(fitted$train_matrix)
  )

  tibble::tibble(
    distance = distance,
    method = method,
    fit = list(fitted),
    cluster = list(assignment),
    average_silhouette = mean(silhouette[, "sil_width"])
  )
}

set.seed(2026)

clustering_results <- tidyr::expand_grid(
  distance = names(distance_presets),
  method = c("PAM", "Spectral")
) |>
  purrr::pmap_dfr(fit_clustering_pipeline) |>
  dplyr::mutate(
    distance_label = unname(distance_labels[distance])
  )

clustering_summary <- clustering_results |>
  dplyr::select(
    distance,
    distance_label,
    method,
    average_silhouette,
    cluster
  ) |>
  tidyr::pivot_wider(
    names_from = method,
    values_from = c(average_silhouette, cluster),
    names_sep = "_"
  ) |>
  dplyr::mutate(
    method_ari = purrr::map2_dbl(
      cluster_PAM,
      cluster_Spectral,
      aricode::ARI
    )
  ) |>
  dplyr::arrange(match(distance, names(distance_presets)))

knitr::kable(
  clustering_summary |>
    dplyr::transmute(
      Distance = distance_label,
      `PAM silhouette` = average_silhouette_PAM,
      `Spectral silhouette` = average_silhouette_Spectral,
      `PAM--spectral ARI` = method_ari
    ),
  digits = 3
)

wdi_benchmark_specs <- all_dist_method_specs(
  mode = "presets_only",
  preset = unname(distance_presets)
) |>
  dplyr::mutate(
    label = unname(distance_labels[preset])
  ) |>
  dplyr::arrange(match(preset, unname(distance_presets)))

set.seed(2026)

wdi_distance_benchmark <- benchmark_mdist(
  x = wdi_clustering,
  specs = wdi_benchmark_specs,
  cluster_k = 4,
  cluster_methods = c("pam", "spectral"),
  spectral_nstart = 50
)

ggplot2::autoplot(
  wdi_distance_benchmark,
  metric = "ari",
  digits = 2
) +
  ggplot2::labs(fill = "ARI") +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold")
  )

### Supervised income classification

## wdi_classification, wdi_split, wdi_train, wdi_folds and neighbor_grid
## were created above, before the refitting comparison.

## The recipes contain an outcome, so u_dep and u_mix use the income
## labels of the current analysis fold; the other three presets are
## response-independent.

classification_recipes <- purrr::map(
  distance_presets,
  ~ recipes::recipe(income ~ ., data = wdi_train) |>
    step_mdist(
      recipes::all_predictors(),
      preset = .x,
      output = "distance_to_training"
    )
)

knn_spec <- nearest_neighbor_dist(
  mode = "classification",
  neighbors = tune::tune()
) |>
  parsnip::set_engine("manydist")

classification_workflows <- workflowsets::workflow_set(
  preproc = classification_recipes,
  models = list(knn = knn_spec)
)

set.seed(2026)

classification_results <- workflowsets::workflow_map(
  classification_workflows,
  fn = "tune_grid",
  seed = 2026,
  resamples = wdi_folds,
  grid = neighbor_grid,
  metrics = yardstick::metric_set(
    yardstick::bal_accuracy,
    yardstick::accuracy
  ),
  control = tune::control_grid(save_pred = TRUE),
  verbose = FALSE
)

classification_ranking <- workflowsets::rank_results(
  classification_results,
  rank_metric = "bal_accuracy",
  select_best = TRUE
)

classification_metrics <- purrr::map_dfr(
  classification_results$wflow_id,
  function(wflow_id) {
    workflowsets::extract_workflow_set_result(
      classification_results,
      wflow_id
    ) |>
      tune::collect_metrics() |>
      dplyr::mutate(wflow_id = wflow_id)
  }
)

classification_summary <- classification_ranking |>
  dplyr::filter(.metric == "bal_accuracy") |>
  dplyr::select(
    wflow_id,
    .config,
    mean,
    std_err,
    rank
  ) |>
  dplyr::left_join(
    classification_metrics |>
      dplyr::filter(.metric == "bal_accuracy") |>
      dplyr::select(
        wflow_id,
        .config,
        neighbors
      ),
    by = c("wflow_id", ".config")
  ) |>
  dplyr::mutate(
    distance = sub("_knn$", "", wflow_id),
    distance_label = unname(distance_labels[distance])
  ) |>
  dplyr::arrange(rank)

knitr::kable(
  classification_summary |>
    dplyr::transmute(
      Distance = distance_label,
      Neighbours = neighbors,
      `Balanced accuracy` = mean,
      `Standard error` = std_err,
      Rank = rank
    ),
  digits = 3
)

## Check. The tuning grid and the refitting comparison estimate the same
## quantity, so the five-neighbour entries below should match the
## `pipeline` column of leakage_comparison.

classification_metrics |>
  dplyr::filter(
    .metric == "bal_accuracy",
    neighbors == 5L,
    wflow_id %in% paste0(leakage_presets, "_knn")
  ) |>
  dplyr::select(wflow_id, neighbors, mean)

classification_metrics |>
  dplyr::filter(.metric == "bal_accuracy") |>
  dplyr::mutate(
    distance = sub("_knn$", "", wflow_id),
    distance_label = unname(distance_labels[distance])
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = neighbors,
      y = mean,
      color = distance_label
    )
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_x_continuous(
    breaks = neighbor_grid$neighbors
  ) +
  ggplot2::scale_color_brewer(
    palette = "Dark2",
    name = "Distance"
  ) +
  ggplot2::labs(
    x = "Number of neighbours",
    y = "Cross-validated balanced accuracy"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  ggplot2::theme(legend.position = "top")

selected_workflow_id <- classification_ranking |>
  dplyr::filter(.metric == "bal_accuracy") |>
  dplyr::slice_min(
    rank,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::pull(wflow_id)

selected_results <- workflowsets::extract_workflow_set_result(
  classification_results,
  selected_workflow_id
)

best_neighbors <- tune::select_best(
  selected_results,
  metric = "bal_accuracy"
)

final_wdi_workflow <- workflowsets::extract_workflow(
  classification_results,
  selected_workflow_id
) |>
  tune::finalize_workflow(best_neighbors)

final_wdi_results <- tune::last_fit(
  final_wdi_workflow,
  split = wdi_split,
  metrics = yardstick::metric_set(
    yardstick::bal_accuracy,
    yardstick::accuracy
  )
)

final_wdi_metrics <- tune::collect_metrics(final_wdi_results)

final_wdi_metrics

selected_distance_label <- classification_summary |>
  dplyr::filter(wflow_id == selected_workflow_id) |>
  dplyr::pull(distance_label)

heldout_balanced_accuracy <- final_wdi_metrics |>
  dplyr::filter(.metric == "bal_accuracy") |>
  dplyr::pull(.estimate)

heldout_accuracy <- final_wdi_metrics |>
  dplyr::filter(.metric == "accuracy") |>
  dplyr::pull(.estimate)

cat("Selected distance:", selected_distance_label, "\n")
cat("Selected neighbours:", best_neighbors$neighbors, "\n")
cat("Held-out balanced accuracy:",
    sprintf("%.3f", heldout_balanced_accuracy), "\n")
cat("Held-out accuracy:",
    sprintf("%.3f", heldout_accuracy), "\n")

wdi_test_predictions <- tune::collect_predictions(
  final_wdi_results
)

wdi_confusion <- yardstick::conf_mat(
  wdi_test_predictions,
  truth = income,
  estimate = .pred_class
)

income_abbreviations <- c(
  "High income" = "High",
  "Low income" = "Low",
  "Lower middle income" = "Lower middle",
  "Upper middle income" = "Upper middle"
)

confusion_table <- wdi_confusion$table |>
  as.data.frame() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    Prediction = dplyr::recode(
      as.character(Prediction),
      !!!income_abbreviations
    ),
    Truth = dplyr::recode(
      as.character(Truth),
      !!!income_abbreviations
    )
  ) |>
  tidyr::pivot_wider(
    names_from = Truth,
    values_from = Freq
  )

knitr::kable(confusion_table)

###################################################
### Computational details
###################################################

sessionInfo()
