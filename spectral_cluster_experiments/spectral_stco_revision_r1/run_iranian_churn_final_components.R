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
      exclusions = c("age_group", "charge_amount")
    )
  )
}

prepared <- r1rd_prepare_data(r1rd_load_iranian_churn(), "Iranian_churn")
data <- prepared$x
truth <- prepared$truth
seed <- analysis_parameters$execution$seed_base + 8000L
screen_file <- file.path(
  raw_data_directory, "Iranian_churn_n3150_screen_499.rds"
)
if (!file.exists(screen_file)) stop("Run the full gated screen first.")
screen <- readRDS(screen_file)
if (screen$n != nrow(data)) stop("The saved gated screen has another sample size.")

output_file <- file.path(raw_data_directory, "Iranian_churn_final_components.rds")
output <- list(
  dataset = "Iranian_churn",
  n = nrow(data),
  class_counts = table(truth),
  roles = list(
    numerical = names(data)[vapply(data, is.numeric, logical(1))],
    categorical = names(data)[vapply(data, is.factor, logical(1))],
    excluded_ordinal = c("age_group", "charge_amount")
  ),
  gated_screen = screen,
  classification = NULL,
  benchmarks = data.frame(
    method = c("D_m", "D_mint"),
    family = "spectral",
    ARI = c(screen$ARI_D_m, screen$ARI_D_mint),
    elapsed_seconds = NA_real_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
)
saveRDS(output, output_file)

cat("Running repeated cross-validated classification diagnostic...\n")
flush.console()
classification_specification$interaction_terms <- NULL
classification <- r1rd_safe_fit(r1rd_run_classification_diagnostic(
  data, truth, seed + classification_specification$seed_offset
))
output$classification <- classification
saveRDS(output, output_file)
if (is.null(classification$value)) {
  cat("Classification failed:", classification$error, "\n")
} else {
  print(classification$value$performance)
  print(classification$value$omnibus)
}

append_benchmark <- function(method, family, expression) {
  cat("Running", method, "...\n")
  flush.console()
  fit <- r1rd_safe_fit(expression)
  row <- data.frame(
    method = method,
    family = family,
    ARI = if (is.null(fit$value)) NA_real_ else fit$value$ARI,
    elapsed_seconds = fit$seconds,
    error = fit$error,
    stringsAsFactors = FALSE
  )
  output$benchmarks <<- rbind(output$benchmarks, row)
  saveRDS(output, output_file)
  print(row)
  invisible(gc())
}

append_benchmark("Gower", "spectral", {
  distance <- as.matrix(manydist::mdist(data, preset = "gower")$distance)
  r1rd_cluster_distance(distance, truth, seed)
})

append_benchmark("Modified Gower", "spectral", {
  distance <- as.matrix(manydist::mdist(data, preset = "mod_gower")$distance)
  r1rd_cluster_distance(distance, truth, seed)
})

append_benchmark("Euclidean one-hot", "spectral", {
  distance <- as.matrix(manydist::mdist(data, preset = "euclidean")$distance)
  r1rd_cluster_distance(distance, truth, seed)
})

append_benchmark("Ahmad--Dey", "spectral", {
  distance <- r1rd_ahmad_dey_distance(
    data, bins = analysis_parameters$distance$ahmad_dey_bins
  )
  r1rd_cluster_distance(distance, truth, seed)
})

append_benchmark("k-prototypes", "algorithm", {
  r1rd_run_kprototypes(data, truth, seed)
})

append_benchmark("KAMILA", "algorithm", {
  r1rd_run_kamila(data, truth, seed)
})

cat("Final benchmark table:\n")
print(output$benchmarks)
cat("Saved:", output_file, "\n")
