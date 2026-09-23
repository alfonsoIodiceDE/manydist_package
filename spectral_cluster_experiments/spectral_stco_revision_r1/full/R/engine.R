full_validate_distance <- function(x, label = "Distance") {
  x <- as.matrix(x)
  if (nrow(x) != ncol(x) || any(!is.finite(x))) {
    stop(label, " must be a finite square matrix.", call. = FALSE)
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-9))) {
    stop(label, " must be symmetric.", call. = FALSE)
  }
  if (any(x < -sqrt(.Machine$double.eps))) {
    stop(label, " contains negative values.", call. = FALSE)
  }
  x <- (x + t(x)) / 2
  x[x < 0] <- 0
  diag(x) <- 0
  x
}

full_mean_distinct <- function(x) {
  x <- as.matrix(x)
  if (nrow(x) < 2L) return(NA_real_)
  mean(x[upper.tri(x)])
}

full_affinity <- function(distance, name, specification) {
  distance <- full_validate_distance(distance)
  n <- nrow(distance)
  if (identical(name, "gaussian")) {
    positive <- distance[upper.tri(distance) & distance > 0]
    sigma <- stats::median(positive)
    if (!length(positive) || !is.finite(sigma) || sigma <= 0) {
      stop("Gaussian affinity has no positive finite scale.", call. = FALSE)
    }
    affinity <- exp(-(distance^2) / (2 * sigma^2))
    diag(affinity) <- 0
    return(list(matrix = affinity, scale = sigma, neighbor_rank = NA_integer_))
  }

  if (identical(name, "selftune")) {
    neighbor_rank <- as.integer(specification$neighbor_rank)
    local_scale <- apply(distance, 1L, function(row) {
      positive <- sort(row[is.finite(row) & row > 0])
      if (!length(positive)) return(NA_real_)
      positive[min(neighbor_rank, length(positive))]
    })
    if (any(!is.finite(local_scale)) || any(local_scale <= 0)) {
      stop("Self-tuning affinity has an invalid local scale.", call. = FALSE)
    }
    denominator <- outer(local_scale, local_scale)
    affinity <- exp(-(distance^2) / denominator)
    diag(affinity) <- 0
    affinity <- (affinity + t(affinity)) / 2
    return(list(
      matrix = affinity,
      scale = stats::median(local_scale),
      neighbor_rank = neighbor_rank
    ))
  }
  stop("Unknown affinity: ", name, call. = FALSE)
}

full_spectral_cluster <- function(affinity, k, seed, kmeans_specification) {
  affinity <- as.matrix(affinity)
  n <- nrow(affinity)
  if (n != ncol(affinity) || k < 2L || k >= n) {
    stop("Invalid affinity dimensions or cluster count.", call. = FALSE)
  }
  degree <- rowSums(affinity)
  if (any(!is.finite(degree)) || any(degree <= 0)) {
    stop("Affinity graph contains a non-positive degree.", call. = FALSE)
  }
  inverse_root_degree <- 1 / sqrt(degree)
  normalized_adjacency <- affinity * outer(inverse_root_degree, inverse_root_degree)
  laplacian <- diag(n) - normalized_adjacency
  decomposition <- eigen(laplacian, symmetric = TRUE)
  selected <- order(decomposition$values)[seq_len(k)]
  embedding <- decomposition$vectors[, selected, drop = FALSE]
  row_norm <- sqrt(rowSums(embedding^2))
  if (any(!is.finite(row_norm)) || any(row_norm <= 0)) {
    stop("Spectral embedding contains an invalid row norm.", call. = FALSE)
  }
  embedding <- embedding / row_norm

  set.seed(seed)
  fit <- stats::kmeans(
    embedding,
    centers = k,
    nstart = as.integer(kmeans_specification$nstart),
    iter.max = as.integer(kmeans_specification$iter_max),
    algorithm = as.character(kmeans_specification$algorithm)
  )
  list(
    cluster = fit$cluster,
    objective = fit$tot.withinss,
    iterations = fit$iter,
    ifault = fit$ifault %||% NA_integer_
  )
}

full_evaluation_settings <- function(task, design) {
  default_prop <- as.numeric(design$distance$defaults$prop_nn)
  default_gamma <- as.numeric(design$distance$defaults$gamma)
  base_prop <- as.numeric(task$parameters$prop_nn %||% default_prop)
  base_gamma <- as.numeric(task$parameters$gamma %||% default_gamma)
  base_numeric_weight <- as.numeric(task$parameters$numeric_block_weight %||% 1)

  within <- task$within_task_grid
  if (!ncol(within)) {
    return(data.frame(
      setting_id = 1L,
      prop_nn = base_prop,
      gamma = base_gamma,
      numeric_block_weight = base_numeric_weight,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    setting_id = seq_len(nrow(within)),
    prop_nn = if ("prop_nn" %in% names(within)) {
      as.numeric(within$prop_nn)
    } else {
      rep(base_prop, nrow(within))
    },
    gamma = if ("gamma" %in% names(within)) {
      as.numeric(within$gamma)
    } else {
      rep(base_gamma, nrow(within))
    },
    numeric_block_weight = if ("numeric_block_weight" %in% names(within)) {
      as.numeric(within$numeric_block_weight)
    } else {
      rep(base_numeric_weight, nrow(within))
    },
    stringsAsFactors = FALSE
  )
}

full_enabled_methods <- function(design, group) {
  selected <- design$methods[vapply(design$methods, function(specification) {
    isTRUE(specification$enabled) && group %in% unlist(specification$groups)
  }, logical(1))]
  selected
}

full_enabled_affinities <- function(design) {
  design$spectral$affinities[vapply(
    design$spectral$affinities,
    function(specification) isTRUE(specification$enabled),
    logical(1)
  )]
}

full_component_level_mean <- function(x) {
  if (nrow(x) < 2L) return(NA_real_)
  mean(x[upper.tri(x)])
}

full_prop_key <- function(prop) {
  format(prop, digits = 15, scientific = FALSE, trim = TRUE)
}

full_prepare_r1_cache <- function(data, settings, design, repo_root, gate_seed) {
  default_prop <- as.numeric(design$distance$defaults$prop_nn)
  props <- unique(if (nrow(settings)) settings$prop_nn else default_prop)
  gate_specification <- design$distance$categorical$interaction_gate
  fits_by_prop <- list()
  seconds_by_prop <- numeric()
  provenance <- scmix_manydist_provenance(
    manydist_repo = repo_root,
    strict = isTRUE(design$manydist$strict_provenance)
  )

  for (prop in props) {
    key <- full_prop_key(prop)
    start <- proc.time()[[3L]]
    fits_by_prop[[key]] <- scmix_gated_active_distance(
      data = data,
      gamma = 1,
      prop_nn = prop,
      score = as.character(design$distance$defaults$score),
      decision = as.character(design$distance$defaults$decision),
      gate_permutations = as.integer(gate_specification$permutations),
      gate_alpha = as.numeric(gate_specification$alpha),
      gate_seed = as.integer(gate_seed)
    )
    seconds_by_prop[[key]] <- proc.time()[[3L]] - start
  }

  base_fit <- fits_by_prop[[full_prop_key(props[[1L]])]]
  list(
    fit = base_fit,
    fits_by_prop = fits_by_prop,
    seconds_by_prop = seconds_by_prop,
    default_gamma = as.numeric(design$distance$defaults$gamma),
    provenance = provenance
  )
}

full_r1_distances_for_setting <- function(cache, data, setting) {
  prop <- setting$prop_nn
  key <- full_prop_key(prop)
  fit <- cache$fits_by_prop[[key]]
  if (is.null(fit)) stop("R1 fit cache miss for prop_nn.", call. = FALSE)
  gamma <- as.numeric(setting$gamma %||% cache$default_gamma)
  if (!is.finite(gamma) || gamma < 0 || gamma > 1) {
    stop("Sensitivity `gamma` must be in [0, 1].", call. = FALSE)
  }
  numeric_weight <- as.numeric(setting$numeric_block_weight)
  if (!is.finite(numeric_weight) || numeric_weight <= 0) {
    stop("Sensitivity `numeric_block_weight` must be positive.", call. = FALSE)
  }
  numeric_distance <- fit$components$numeric
  categorical_distance <- fit$components$categorical
  active_average <- fit$components$interaction_active_average
  weighted_numeric <- numeric_weight * numeric_distance
  baseline <- full_validate_distance(
    weighted_numeric + categorical_distance,
    "Pure package no-interaction baseline"
  )
  interaction_scale <- gamma * full_mean_distinct(weighted_numeric)
  added_interaction <- interaction_scale * active_average
  interaction_distance <- if (!fit$diagnostics$active_count || gamma == 0) {
    baseline
  } else {
    full_validate_distance(
      baseline + added_interaction,
      "Gated active-average interaction distance"
    )
  }
  if (!fit$diagnostics$active_count && !identical(interaction_distance, baseline)) {
    stop("No active gate must recover the package baseline exactly.", call. = FALSE)
  }

  variables <- names(fit$components$interaction_level_delta)
  diagnostic_rows <- lapply(seq_along(variables), function(j) {
    variable <- variables[[j]]
    interaction_raw <- fit$components$interaction_observation[[variable]]
    data.frame(
      variable = variable,
      prop_nn = prop,
      gamma = gamma,
      numeric_block_weight = numeric_weight,
      tvd_level_mean = NA_real_,
      interaction_level_mean = full_component_level_mean(
        fit$components$interaction_level_delta[[variable]]
      ),
      tvd_observation_mean = NA_real_,
      interaction_observation_mean_raw = full_mean_distinct(
        interaction_raw
      ),
      interaction_observation_mean_gated = if (isTRUE(fit$gate$detected[[variable]])) {
        full_mean_distinct(interaction_raw) / fit$diagnostics$active_count
      } else {
        0
      },
      interaction_gate_statistic = fit$gate$statistic[[variable]],
      interaction_gate_p_value = fit$gate$p_value[[variable]],
      interaction_detected = isTRUE(fit$gate$detected[[variable]]),
      stringsAsFactors = FALSE
    )
  })

  list(
    distances = list(
      Gated_active_average_SC = interaction_distance,
      Package_no_interaction_SC = baseline
    ),
    prop_nn = prop,
    gamma = gamma,
    numeric_block_weight = numeric_weight,
    diagnostics = do.call(rbind, diagnostic_rows),
    active_count = fit$diagnostics$active_count,
    active_names = fit$diagnostics$active_names,
    interaction_scale = interaction_scale,
    added_interaction_mean = full_mean_distinct(added_interaction),
    distance_seconds = cache$seconds_by_prop[[key]]
  )
}

full_manydist_distance <- function(data, engine, prop_nn, design) {
  if (identical(engine, "manydist_gower_distance")) {
    return(as.matrix(manydist::mdist(data, preset = "gower")$distance))
  }
  if (identical(engine, "manydist_euclidean_distance")) {
    return(as.matrix(manydist::mdist(data, preset = "euclidean")$distance))
  }
  if (identical(engine, "manydist_modified_gower_distance")) {
    return(as.matrix(manydist::mdist(data, preset = "mod_gower")$distance))
  }
  if (identical(engine, "manydist_dkss_distance")) {
    return(as.matrix(manydist::mdist(data, preset = "dkss")$distance))
  }
  if (identical(engine, "legacy_u_dep_interaction_distance")) {
    fit <- manydist::mdist(
      data,
      preset = "custom",
      method_cat = "tvd",
      method_num = "pc_scores",
      commensurable = FALSE,
      interaction = TRUE,
      prop_nn = prop_nn,
      score = as.character(design$distance$defaults$score),
      decision = as.character(design$distance$defaults$decision)
    )
    return(as.matrix(fit$distance))
  }
  stop("Unknown manydist engine: ", engine, call. = FALSE)
}

full_score_distance <- function(
    distance,
    truth,
    affinity_name,
    affinity_specification,
    cluster_seed,
    design
) {
  affinity <- full_affinity(distance, affinity_name, affinity_specification)
  fit <- full_spectral_cluster(
    affinity$matrix,
    k = nlevels(factor(truth)),
    seed = cluster_seed,
    kmeans_specification = design$spectral$kmeans
  )
  list(
    ARI = mclust::adjustedRandIndex(fit$cluster, truth),
    objective = fit$objective,
    iterations = fit$iterations,
    ifault = fit$ifault,
    affinity_scale = affinity$scale,
    affinity_neighbor_rank = affinity$neighbor_rank
  )
}

full_result_row <- function(
    task, setting, method_name, method_label, affinity, n, k,
    prop_nn, gamma, numeric_block_weight, distance_seconds, cluster_seconds,
    score = NULL, error = NA_character_
) {
  data.frame(
    task_id = task$task_id,
    scenario = task$scenario,
    condition_id = task$condition_id,
    replicate = task$replicate,
    data_seed = task$data_seed,
    setting_id = setting$setting_id,
    n = n,
    k = k,
    prop_nn = prop_nn,
    gamma = gamma,
    numeric_block_weight = numeric_block_weight,
    method = method_name,
    method_label = method_label,
    affinity = affinity,
    ARI = if (is.null(score)) NA_real_ else score$ARI,
    objective = if (is.null(score)) NA_real_ else score$objective,
    iterations = if (is.null(score)) NA_integer_ else score$iterations,
    ifault = if (is.null(score)) NA_integer_ else score$ifault,
    affinity_scale = if (is.null(score)) NA_real_ else score$affinity_scale,
    affinity_neighbor_rank = if (is.null(score)) NA_integer_ else score$affinity_neighbor_rank,
    distance_seconds = distance_seconds,
    clustering_seconds = cluster_seconds,
    error = error,
    stringsAsFactors = FALSE
  )
}

full_run_kprototypes <- function(data, truth, seed, design) {
  if (!requireNamespace("clustMixType", quietly = TRUE)) {
    stop("The enabled KPrototypes method requires `clustMixType`.", call. = FALSE)
  }
  set.seed(seed)
  fit <- clustMixType::kproto(
    data,
    k = nlevels(factor(truth)),
    lambda = NULL,
    type = "huang",
    iter.max = as.integer(design$spectral$kmeans$iter_max),
    nstart = as.integer(design$spectral$kmeans$nstart),
    na.rm = "no",
    keep.data = FALSE,
    verbose = FALSE
  )
  list(
    ARI = mclust::adjustedRandIndex(fit$cluster, truth),
    objective = fit$tot.withinss %||% fit$withinss %||% NA_real_,
    iterations = fit$iter %||% NA_integer_,
    ifault = NA_integer_,
    affinity_scale = NA_real_,
    affinity_neighbor_rank = NA_integer_,
    lambda = fit$lambda %||% NA_real_
  )
}

full_evaluate_task <- function(task, generated, design, repo_root) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("The experiment runner requires `mclust` for ARI.", call. = FALSE)
  }
  data <- generated$x
  truth <- droplevels(factor(generated$truth))
  settings <- full_evaluation_settings(task, design)
  methods <- full_enabled_methods(design, task$method_group)
  affinities <- full_enabled_affinities(design)
  method_names <- names(methods)
  result_rows <- list()
  diagnostic_rows <- list()
  result_index <- 0L
  diagnostic_index <- 0L

  needs_r1 <- any(vapply(methods, function(x) {
    x$engine %in% c("r1_interaction_distance", "r1_no_interaction_distance")
  }, logical(1)))
  r1_error <- NA_character_
  r1_cache <- NULL
  if (needs_r1) {
    prepared <- tryCatch(
      list(
        value = full_prepare_r1_cache(
          data,
          settings,
          design,
          repo_root,
          gate_seed = task$data_seed +
            as.integer(design$distance$categorical$interaction_gate$seed_offset)
        ),
        error = NA_character_
      ),
      error = function(e) list(value = NULL, error = conditionMessage(e))
    )
    r1_cache <- prepared$value
    r1_error <- prepared$error
  }

  # Core comparator distances do not vary over a within-task sensitivity grid.
  comparator_distances <- list()
  comparator_seconds <- numeric()
  comparator_errors <- list()
  for (method_name in method_names) {
    engine <- methods[[method_name]]$engine
    if (!engine %in% c(
      "legacy_u_dep_interaction_distance",
      "manydist_gower_distance",
      "manydist_euclidean_distance",
      "manydist_modified_gower_distance",
      "manydist_dkss_distance"
    )) next

    start <- proc.time()[[3L]]
    fitted <- tryCatch(
      list(
        distance = full_validate_distance(
          full_manydist_distance(data, engine, settings$prop_nn[[1L]], design),
          method_name
        ),
        error = NA_character_
      ),
      error = function(e) list(distance = NULL, error = conditionMessage(e))
    )
    comparator_seconds[[method_name]] <- proc.time()[[3L]] - start
    comparator_distances[[method_name]] <- fitted$distance
    comparator_errors[[method_name]] <- fitted$error
  }

  for (setting_row in seq_len(nrow(settings))) {
    setting <- settings[setting_row, , drop = FALSE]
    r1_setting_error <- r1_error
    r1_setting <- NULL
    if (needs_r1 && !is.null(r1_cache)) {
      prepared_setting <- tryCatch(
        list(
          value = full_r1_distances_for_setting(r1_cache, data, setting),
          error = NA_character_
        ),
        error = function(e) list(value = NULL, error = conditionMessage(e))
      )
      r1_setting <- prepared_setting$value
      r1_setting_error <- prepared_setting$error
    }
    if (!is.null(r1_setting)) {
      diagnostics <- r1_setting$diagnostics
      diagnostics$task_id <- task$task_id
      diagnostics$scenario <- task$scenario
      diagnostics$replicate <- task$replicate
      diagnostics$setting_id <- setting$setting_id
      diagnostic_index <- diagnostic_index + 1L
      diagnostic_rows[[diagnostic_index]] <- diagnostics
    }

    for (method_name in method_names) {
      specification <- methods[[method_name]]
      engine <- specification$engine

      if (identical(engine, "clustMixType_kprototypes")) {
        start <- proc.time()[[3L]]
        fitted <- tryCatch(
          list(
            score = full_run_kprototypes(
              data,
              truth,
              task$data_seed + as.integer(design$execution$cluster_seed_offsets$kprototypes),
              design
            ),
            error = NA_character_
          ),
          error = function(e) list(score = NULL, error = conditionMessage(e))
        )
        elapsed <- proc.time()[[3L]] - start
        result_index <- result_index + 1L
        result_rows[[result_index]] <- full_result_row(
          task, setting, method_name, specification$label,
          affinity = "not_applicable", n = nrow(data), k = nlevels(truth),
          prop_nn = r1_setting$prop_nn %||% setting$prop_nn,
          gamma = r1_setting$gamma %||% setting$gamma,
          numeric_block_weight = r1_setting$numeric_block_weight %||% setting$numeric_block_weight,
          distance_seconds = NA_real_, cluster_seconds = elapsed,
          score = fitted$score, error = fitted$error
        )
        next
      }

      if (identical(engine, "r1_interaction_distance")) {
        distance <- if (is.null(r1_setting)) NULL else {
          r1_setting$distances$Gated_active_average_SC
        }
        distance_seconds <- if (is.null(r1_setting)) NA_real_ else r1_setting$distance_seconds
        distance_error <- r1_setting_error
      } else if (identical(engine, "r1_no_interaction_distance")) {
        distance <- if (is.null(r1_setting)) NULL else {
          r1_setting$distances$Package_no_interaction_SC
        }
        distance_seconds <- if (is.null(r1_setting)) NA_real_ else r1_setting$distance_seconds
        distance_error <- r1_setting_error
      } else {
        distance <- comparator_distances[[method_name]]
        distance_seconds <- comparator_seconds[[method_name]] %||% NA_real_
        distance_error <- comparator_errors[[method_name]] %||% "Distance was not constructed."
      }

      for (affinity_name in names(affinities)) {
        if (is.null(distance)) {
          result_index <- result_index + 1L
          result_rows[[result_index]] <- full_result_row(
            task, setting, method_name, specification$label,
            affinity = affinity_name, n = nrow(data), k = nlevels(truth),
            prop_nn = r1_setting$prop_nn %||% setting$prop_nn,
            gamma = r1_setting$gamma %||% setting$gamma,
            numeric_block_weight = r1_setting$numeric_block_weight %||% setting$numeric_block_weight,
            distance_seconds = distance_seconds,
            cluster_seconds = NA_real_, error = distance_error
          )
          next
        }

        cluster_seed <- task$data_seed +
          as.integer(design$execution$cluster_seed_offsets[[affinity_name]])
        start <- proc.time()[[3L]]
        scored <- tryCatch(
          list(
            score = full_score_distance(
              distance,
              truth,
              affinity_name,
              affinities[[affinity_name]],
              cluster_seed,
              design
            ),
            error = NA_character_
          ),
          error = function(e) list(score = NULL, error = conditionMessage(e))
        )
        cluster_seconds <- proc.time()[[3L]] - start
        result_index <- result_index + 1L
        result_rows[[result_index]] <- full_result_row(
          task, setting, method_name, specification$label,
          affinity = affinity_name, n = nrow(data), k = nlevels(truth),
          prop_nn = r1_setting$prop_nn %||% setting$prop_nn,
          gamma = r1_setting$gamma %||% setting$gamma,
          numeric_block_weight = r1_setting$numeric_block_weight %||% setting$numeric_block_weight,
          distance_seconds = distance_seconds,
          cluster_seconds = cluster_seconds,
          score = scored$score, error = scored$error
        )
      }
    }
  }

  list(
    results = do.call(rbind, result_rows),
    diagnostics = if (length(diagnostic_rows)) do.call(rbind, diagnostic_rows) else data.frame(),
    r1_provenance = if (is.null(r1_cache)) NULL else r1_cache$provenance,
    r1_error = r1_error
  )
}
