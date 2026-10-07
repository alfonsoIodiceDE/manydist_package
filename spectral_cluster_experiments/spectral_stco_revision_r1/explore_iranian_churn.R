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

r1rd_load_iranian_churn <- function() {
  path <- file.path(
    raw_data_directory, "iranian_churn_source", "Customer Churn.csv"
  )
  raw <- utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE
  )
  names(raw) <- r1rd_normalize_name(names(raw))
  expected <- c(
    "call_failure", "complains", "subscription_length", "charge_amount",
    "seconds_of_use", "frequency_of_use", "frequency_of_sms",
    "distinct_called_numbers", "age_group", "tariff_plan", "status", "age",
    "customer_value", "churn"
  )
  if (!identical(names(raw), expected)) {
    stop("Iranian Churn source does not contain the expected columns.")
  }

  numeric_names <- c(
    "call_failure", "subscription_length", "seconds_of_use",
    "frequency_of_use", "frequency_of_sms", "distinct_called_numbers",
    "age", "customer_value"
  )
  categorical_names <- c("complains", "tariff_plan", "status")
  raw[numeric_names] <- lapply(raw[numeric_names], as.numeric)
  raw[categorical_names] <- lapply(raw[categorical_names], factor)

  list(
    x = raw[c(numeric_names, categorical_names)],
    truth = factor(raw$churn, levels = c(0, 1), labels = c("retained", "churned")),
    provenance = list(
      source = "UCI Iranian Churn",
      doi = "10.24432/C5JW3Z",
      file = path,
      exclusions = paste(
        "age_group excluded as a coarsened duplicate of age;",
        "ordinal charge_amount excluded because ordinal dissimilarities",
        "are outside the present nominal--numerical scope"
      )
    )
  )
}

prepared <- r1rd_prepare_data(r1rd_load_iranian_churn(), "Iranian_churn")
data <- prepared$x
truth <- prepared$truth
seed <- analysis_parameters$execution$seed_base + 8000L
arguments <- commandArgs(trailingOnly = TRUE)
gate_permutations <- as.integer(arguments[1L] %||% 99L)
if (is.na(gate_permutations)) gate_permutations <- 99L
sample_size <- as.integer(arguments[2L] %||% nrow(data))
if (is.na(sample_size) || sample_size < 1L) sample_size <- nrow(data)

full_n <- nrow(data)
full_class_counts <- table(truth)
if (sample_size < full_n) {
  proportions <- as.numeric(full_class_counts) / sum(full_class_counts)
  allocation <- floor(sample_size * proportions)
  remainder <- sample_size - sum(allocation)
  if (remainder > 0L) {
    priorities <- order(sample_size * proportions - allocation, decreasing = TRUE)
    allocation[priorities[seq_len(remainder)]] <-
      allocation[priorities[seq_len(remainder)]] + 1L
  }
  sampled_rows <- .scmix_with_seed(seed, unlist(Map(
    function(level, size) sample(which(truth == level), size = size),
    levels(truth), allocation
  )))
  sampled_rows <- sort(sampled_rows)
  data <- data[sampled_rows, , drop = FALSE]
  truth <- droplevels(truth[sampled_rows])
}

cat(
  "Iranian Churn:", nrow(data), "observations;",
  sum(vapply(data, is.numeric, logical(1))), "numerical and",
  sum(vapply(data, is.factor, logical(1))), "categorical predictors.\n"
)
cat("Computing gated distances with", gate_permutations, "permutations...\n")
flush.console()

started <- Sys.time()
fit <- scmix_gated_active_distance(
  data = data,
  gamma = analysis_parameters$distance$gamma,
  prop_nn = analysis_parameters$distance$prop_nn,
  score = analysis_parameters$distance$score,
  decision = analysis_parameters$distance$decision,
  gate_permutations = gate_permutations,
  gate_alpha = analysis_parameters$distance$gate_alpha,
  gate_seed = seed + analysis_parameters$distance$gate_seed_offset
)
cat("Gate complete; computing the two spectral partitions...\n")
flush.console()
dm <- r1rd_cluster_distance(fit$distance_no_interaction, truth, seed)
dmint <- r1rd_cluster_distance(fit$distance, truth, seed)

diagnostics <- list(
  full_n = full_n,
  full_class_counts = full_class_counts,
  n = nrow(data),
  class_counts = table(truth),
  p_numeric = sum(vapply(data, is.numeric, logical(1))),
  q_categorical = sum(vapply(data, is.factor, logical(1))),
  excluded_variables = c("age_group", "charge_amount"),
  gate_permutations = gate_permutations,
  gate = data.frame(
    variable = names(fit$gate$detected),
    detected = unname(fit$gate$detected),
    multiplicity_adjusted_p_value = unname(fit$gate$p_value),
    statistic = unname(fit$gate$statistic)
  ),
  active_variables = fit$diagnostics$active_names,
  ARI_D_m = dm$ARI,
  ARI_D_mint = dmint$ARI,
  partition_ARI = mclust::adjustedRandIndex(dm$cluster, dmint$cluster),
  R_align = if (fit$diagnostics$active_count > 0L) {
    r1rd_alignment_ratio(fit$components$interaction_active_average, dmint$cluster)
  } else {
    NA_real_
  },
  interaction_baseline_spearman = if (fit$diagnostics$active_count > 0L) {
    r1rd_distance_spearman(
      fit$components$interaction_active_average,
      fit$distance_no_interaction
    )
  } else {
    NA_real_
  },
  interaction_numeric_spearman = if (fit$diagnostics$active_count > 0L) {
    r1rd_distance_spearman(
      fit$components$interaction_active_average, fit$components$numeric
    )
  } else {
    NA_real_
  },
  interaction_categorical_spearman = if (fit$diagnostics$active_count > 0L) {
    r1rd_distance_spearman(
      fit$components$interaction_active_average, fit$components$categorical
    )
  } else {
    NA_real_
  },
  neighbourhood_overlap = r1rd_knn_overlap(
    fit$distance_no_interaction, fit$distance,
    analysis_parameters$distance$prop_nn
  ),
  elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))
)

print(diagnostics)
saveRDS(
  diagnostics,
  file.path(
    raw_data_directory,
    paste0(
      "Iranian_churn_n", nrow(data), "_screen_", gate_permutations, ".rds"
    )
  )
)
