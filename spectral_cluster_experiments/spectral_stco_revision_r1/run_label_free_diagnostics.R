qmd <- normalizePath("R1_real_data.qmd", mustWork = TRUE)
lines <- readLines(qmd, warn = FALSE)
lines <- lines[seq_len(match("# Optional execution", lines) - 1L)]

chunks <- list()
inside <- FALSE
current <- character()
for (line in lines) {
  if (!inside && grepl("^```\\{r", line)) {
    inside <- TRUE
    current <- character()
  } else if (inside && grepl("^```[[:space:]]*$", line)) {
    chunks[[length(chunks) + 1L]] <- current
    inside <- FALSE
  } else if (inside && !grepl("^#\\|", line)) {
    current <- c(current, line)
  }
}
for (chunk in chunks) {
  eval(parse(text = paste(chunk, collapse = "\n")), envir = .GlobalEnv)
}

arguments <- commandArgs(trailingOnly = TRUE)
selected_datasets <- strsplit(
  Sys.getenv(
    "R1_DIAGNOSTIC_DATASETS",
    "Star,Palmer_penguins,SAheart,Dermatology"
  ),
  ",",
  fixed = TRUE
)[[1L]]
selected_datasets <- trimws(selected_datasets)
resamples <- as.integer(arguments[1L] %||% 50L)
if (!is.finite(resamples) || resamples < 2L) resamples <- 50L
sample_fraction <- 0.80
bandwidth_multipliers <- c(0.50, 0.75, 1.00, 1.50, 2.00)
gamma_grid <- c(0, 0.25, 0.50, 0.75, 1.00)

output_directory <- file.path(
  revision_directory,
  "results",
  Sys.getenv("R1_DIAGNOSTIC_OUTPUT", "R1_label_free_diagnostics")
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

diagnostic_spectral <- function(distance, k, seed, bandwidth_multiplier = 1) {
  distance <- full_validate_distance(distance)
  positive <- distance[upper.tri(distance) & distance > 0]
  sigma <- stats::median(positive) * bandwidth_multiplier
  if (!length(positive) || !is.finite(sigma) || sigma <= 0) {
    stop("Invalid bandwidth in diagnostic spectral clustering.")
  }
  affinity <- exp(-(distance^2) / (2 * sigma^2))
  diag(affinity) <- 0
  degree <- rowSums(affinity)
  inverse_root_degree <- 1 / sqrt(degree)
  normalized_adjacency <- affinity * outer(
    inverse_root_degree, inverse_root_degree
  )
  laplacian <- diag(nrow(distance)) - normalized_adjacency
  decomposition <- eigen(laplacian, symmetric = TRUE)
  order_index <- order(decomposition$values)
  eigenvalues <- decomposition$values[order_index]
  embedding <- decomposition$vectors[, order_index[seq_len(k)], drop = FALSE]
  embedding <- embedding / sqrt(rowSums(embedding^2))
  set.seed(seed)
  fit <- stats::kmeans(
    embedding,
    centers = k,
    nstart = kmeans_specification$nstart,
    iter.max = kmeans_specification$iter_max,
    algorithm = kmeans_specification$algorithm
  )
  list(
    cluster = fit$cluster,
    eigengap = eigenvalues[[k + 1L]] - eigenvalues[[k]],
    lambda_k = eigenvalues[[k]],
    lambda_k_plus_1 = eigenvalues[[k + 1L]],
    sigma = sigma
  )
}

mean_distinct <- function(distance) {
  mean(distance[upper.tri(distance)])
}

bind_rows <- function(rows) {
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

summarize_stability <- function(values) {
  c(
    mean = mean(values),
    median = stats::median(values),
    q10 = unname(stats::quantile(values, 0.10)),
    q90 = unname(stats::quantile(values, 0.90)),
    sd = stats::sd(values)
  )
}

dataset_outputs <- list()

for (dataset_index in seq_along(selected_datasets)) {
  dataset <- selected_datasets[[dataset_index]]
  message("Label-free diagnostics: ", dataset)
  prepared <- r1rd_prepare_data(r1rd_load_dataset(dataset), dataset)
  data <- prepared$x
  truth <- prepared$truth
  k <- nlevels(truth)
  seed <- analysis_parameters$execution$seed_base +
    match(dataset, dataset_registry$dataset) * 1000L

  gated_fit <- scmix_gated_active_distance(
    data = data,
    gamma = analysis_parameters$distance$gamma,
    prop_nn = analysis_parameters$distance$prop_nn,
    score = analysis_parameters$distance$score,
    decision = analysis_parameters$distance$decision,
    gate_permutations = analysis_parameters$distance$gate_permutations,
    gate_alpha = analysis_parameters$distance$gate_alpha,
    gate_seed = seed + analysis_parameters$distance$gate_seed_offset
  )
  distances <- list(
    D_m = gated_fit$distance_no_interaction,
    D_mint = gated_fit$distance
  )
  full_fits <- lapply(names(distances), function(method) {
    diagnostic_spectral(distances[[method]], k, seed)
  })
  names(full_fits) <- names(distances)

  interaction <- gated_fit$components$interaction_active_average
  added <- gated_fit$diagnostics$added_mean_distinct
  geometry <- data.frame(
    dataset = dataset,
    n = nrow(data),
    K = k,
    active_count = gated_fit$diagnostics$active_count,
    active_variables = paste(gated_fit$diagnostics$active_names, collapse = ", "),
    interaction_baseline_spearman = r1rd_distance_spearman(
      interaction, distances$D_m
    ),
    interaction_numeric_spearman = r1rd_distance_spearman(
      interaction, gated_fit$components$numeric
    ),
    interaction_categorical_spearman = r1rd_distance_spearman(
      interaction, gated_fit$components$categorical
    ),
    neighbourhood_overlap = r1rd_knn_overlap(
      distances$D_m, distances$D_mint,
      analysis_parameters$distance$prop_nn
    ),
    partition_ARI_D_m_vs_D_mint = mclust::adjustedRandIndex(
      full_fits$D_m$cluster, full_fits$D_mint$cluster
    ),
    added_interaction_mean = added,
    added_to_baseline_mean_ratio = added / mean_distinct(distances$D_m),
    stringsAsFactors = FALSE
  )

  primary_spectral <- bind_rows(lapply(names(distances), function(method) {
    fit <- full_fits[[method]]
    data.frame(
      dataset = dataset,
      method = method,
      eigengap = fit$eigengap,
      lambda_k = fit$lambda_k,
      lambda_k_plus_1 = fit$lambda_k_plus_1,
      sigma = fit$sigma,
      stringsAsFactors = FALSE
    )
  }))

  set.seed(seed + 9000L)
  samples <- replicate(
    resamples,
    sort(sample.int(nrow(data), floor(sample_fraction * nrow(data)))),
    simplify = FALSE
  )
  stability_rows <- list()
  row_index <- 0L
  for (method in names(distances)) {
    for (resample_index in seq_along(samples)) {
      index <- samples[[resample_index]]
      subset_fit <- diagnostic_spectral(
        distances[[method]][index, index, drop = FALSE],
        k,
        seed + 10000L + resample_index
      )
      row_index <- row_index + 1L
      stability_rows[[row_index]] <- data.frame(
        dataset = dataset,
        method = method,
        resample = resample_index,
        stability_ARI = mclust::adjustedRandIndex(
          subset_fit$cluster, full_fits[[method]]$cluster[index]
        ),
        eigengap = subset_fit$eigengap,
        stringsAsFactors = FALSE
      )
    }
  }
  stability <- bind_rows(stability_rows)
  stability_summary <- bind_rows(lapply(
    split(stability, stability$method),
    function(observed) {
      ari_summary <- summarize_stability(observed$stability_ARI)
      gap_summary <- summarize_stability(observed$eigengap)
      data.frame(
        dataset = dataset,
        method = observed$method[[1L]],
        stability_mean = ari_summary[["mean"]],
        stability_median = ari_summary[["median"]],
        stability_q10 = ari_summary[["q10"]],
        stability_q90 = ari_summary[["q90"]],
        stability_sd = ari_summary[["sd"]],
        eigengap_resample_mean = gap_summary[["mean"]],
        eigengap_resample_q10 = gap_summary[["q10"]],
        stringsAsFactors = FALSE
      )
    }
  ))

  bandwidth_rows <- list()
  row_index <- 0L
  for (method in names(distances)) {
    for (multiplier in bandwidth_multipliers) {
      fit <- diagnostic_spectral(
        distances[[method]], k,
        seed + as.integer(round(100 * multiplier)),
        bandwidth_multiplier = multiplier
      )
      row_index <- row_index + 1L
      bandwidth_rows[[row_index]] <- data.frame(
        dataset = dataset,
        method = method,
        bandwidth_multiplier = multiplier,
        eigengap = fit$eigengap,
        partition_ARI_vs_primary = mclust::adjustedRandIndex(
          fit$cluster, full_fits[[method]]$cluster
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  bandwidth <- bind_rows(bandwidth_rows)

  gamma_rows <- list()
  for (gamma_index in seq_along(gamma_grid)) {
    gamma <- gamma_grid[[gamma_index]]
    gamma_distance <- full_validate_distance(
      distances$D_m + gamma * mean_distinct(gated_fit$components$numeric) *
        interaction
    )
    fit <- diagnostic_spectral(
      gamma_distance, k, seed + 20000L + gamma_index
    )
    gamma_rows[[gamma_index]] <- data.frame(
      dataset = dataset,
      gamma = gamma,
      eigengap = fit$eigengap,
      partition_ARI_vs_D_m = mclust::adjustedRandIndex(
        fit$cluster, full_fits$D_m$cluster
      ),
      partition_ARI_vs_D_mint = mclust::adjustedRandIndex(
        fit$cluster, full_fits$D_mint$cluster
      ),
      external_ARI = mclust::adjustedRandIndex(fit$cluster, truth),
      stringsAsFactors = FALSE
    )
  }
  gamma_path <- bind_rows(gamma_rows)

  artifact <- readRDS(r1rd_artifact_file(dataset))
  gate_sensitivity <- artifact$neighbour_sensitivity[c(
    "dataset", "prop_nn", "partition_ARI_vs_primary", "active_count",
    "active_variables"
  )]

  external_validation <- data.frame(
    dataset = dataset,
    method = names(full_fits),
    external_ARI = vapply(full_fits, function(fit) {
      mclust::adjustedRandIndex(fit$cluster, truth)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  dataset_outputs[[dataset]] <- list(
    geometry = geometry,
    primary_spectral = primary_spectral,
    stability = stability,
    stability_summary = stability_summary,
    bandwidth = bandwidth,
    gamma_path = gamma_path,
    gate_sensitivity = gate_sensitivity,
    external_validation = external_validation,
    gate = data.frame(
      dataset = dataset,
      variable = names(gated_fit$gate$detected),
      detected = unname(gated_fit$gate$detected),
      adjusted_p = unname(gated_fit$gate$p_value),
      statistic = unname(gated_fit$gate$statistic),
      stringsAsFactors = FALSE
    )
  )
  saveRDS(
    dataset_outputs[[dataset]],
    file.path(output_directory, paste0(dataset, ".rds")),
    version = 3L
  )
}

combined <- list(
  geometry = bind_rows(lapply(dataset_outputs, `[[`, "geometry")),
  primary_spectral = bind_rows(lapply(dataset_outputs, `[[`, "primary_spectral")),
  stability = bind_rows(lapply(dataset_outputs, `[[`, "stability")),
  stability_summary = bind_rows(lapply(dataset_outputs, `[[`, "stability_summary")),
  bandwidth = bind_rows(lapply(dataset_outputs, `[[`, "bandwidth")),
  gamma_path = bind_rows(lapply(dataset_outputs, `[[`, "gamma_path")),
  gate_sensitivity = bind_rows(lapply(dataset_outputs, `[[`, "gate_sensitivity")),
  external_validation = bind_rows(lapply(dataset_outputs, `[[`, "external_validation")),
  gate = bind_rows(lapply(dataset_outputs, `[[`, "gate")),
  specification = list(
    datasets = selected_datasets,
    resamples = resamples,
    sample_fraction = sample_fraction,
    bandwidth_multipliers = bandwidth_multipliers,
    gamma_grid = gamma_grid,
    labels_used_only_for = "fixed K and external validation tables"
  )
)

saveRDS(combined, file.path(output_directory, "combined.rds"), version = 3L)
for (name in setdiff(names(combined), "specification")) {
  utils::write.csv(
    combined[[name]], file.path(output_directory, paste0(name, ".csv")),
    row.names = FALSE
  )
}

print(combined$geometry)
print(combined$primary_spectral)
print(combined$stability_summary)
print(combined$external_validation)
