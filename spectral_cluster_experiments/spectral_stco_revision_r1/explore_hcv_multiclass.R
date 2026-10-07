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

r1rd_load_hcv_multiclass <- function() {
  path <- file.path(raw_data_directory, "HCV.csv")
  raw <- utils::read.csv(
    path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
  )
  names(raw)[1L] <- "patient_id"
  numeric_names <- c(
    "Age", "ALB", "ALP", "ALT", "AST", "BIL", "CHE", "CHOL",
    "CREA", "GGT", "PROT"
  )
  raw[numeric_names] <- lapply(raw[numeric_names], as.numeric)
  diagnosis_labels <- c(
    "blood_donor", "suspected_blood_donor", "hepatitis", "fibrosis",
    "cirrhosis"
  )
  category_codes <- c(
    "0=Blood Donor", "0s=suspect Blood Donor", "1=Hepatitis",
    "2=Fibrosis", "3=Cirrhosis"
  )
  list(
    x = data.frame(raw[numeric_names], Sex = factor(raw$Sex), check.names = FALSE),
    truth = factor(
      diagnosis_labels[match(raw$Category, category_codes)],
      levels = diagnosis_labels
    ),
    provenance = list(
      url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00571/hcvdat0.csv",
      doi = "10.24432/C5D612",
      file = path
    )
  )
}

prepared <- r1rd_prepare_data(r1rd_load_hcv_multiclass(), "HCV")
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
    paste0("HCV_multiclass_screen_", gate_permutations, ".rds")
  )
)
