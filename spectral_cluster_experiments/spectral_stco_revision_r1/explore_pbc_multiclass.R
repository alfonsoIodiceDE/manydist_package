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

r1rd_load_pbc_multiclass <- function() {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The PBC screen requires package `survival`.")
  }
  environment <- new.env(parent = emptyenv())
  utils::data("pbc", package = "survival", envir = environment)
  raw <- as.data.frame(get("pbc", envir = environment, inherits = FALSE))

  numeric_names <- c(
    "age", "bili", "chol", "albumin", "copper", "alk.phos", "ast",
    "trig", "platelet", "protime"
  )
  categorical_names <- c(
    "trt", "sex", "ascites", "hepato", "spiders", "edema", "stage"
  )
  raw$age <- raw$age / 365.25
  raw[numeric_names] <- lapply(raw[numeric_names], as.numeric)
  raw[categorical_names] <- lapply(raw[categorical_names], factor)
  status_labels <- c("censored", "transplant", "death")
  list(
    x = raw[c(numeric_names, categorical_names)],
    truth = factor(status_labels[raw$status + 1L], levels = status_labels),
    provenance = list(
      package = "survival",
      object = "pbc",
      source = "Mayo Clinic primary biliary cirrhosis trial"
    )
  )
}

prepared <- r1rd_prepare_data(r1rd_load_pbc_multiclass(), "PBC_multiclass")
data <- prepared$x
truth <- prepared$truth
seed <- analysis_parameters$execution$seed_base + 7000L
gate_permutations <- as.integer(commandArgs(trailingOnly = TRUE)[1L] %||% 99L)
if (is.na(gate_permutations)) gate_permutations <- 99L

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
dm <- r1rd_cluster_distance(fit$distance_no_interaction, truth, seed)
dmint <- r1rd_cluster_distance(fit$distance, truth, seed)

diagnostics <- list(
  n = nrow(data),
  class_counts = table(truth),
  p_numeric = sum(vapply(data, is.numeric, logical(1))),
  q_categorical = sum(vapply(data, is.factor, logical(1))),
  numeric_missing = prepared$preprocessing$numeric_missing,
  categorical_missing = prepared$preprocessing$categorical_missing,
  gate_permutations = gate_permutations,
  gate = data.frame(
    variable = names(fit$gate$detected),
    detected = unname(fit$gate$detected),
    adjusted_p = unname(fit$gate$p_value),
    statistic = unname(fit$gate$statistic)
  ),
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
    paste0("PBC_multiclass_screen_", gate_permutations, ".rds")
  )
)
