#!/usr/bin/env Rscript

# Small, reproducible pilot for the Statistics & Computing R1 revision.
# This is deliberately not the enriched/full experiment suite.

script_path <- function() {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg)) return(normalizePath(sub("^--file=", "", arg[[1]])))
  normalizePath("spectral_cluster_experiments/revision_r1/pilot/run_pilot.R")
}

pilot_file <- script_path()
revision_dir <- dirname(dirname(pilot_file))
repo_root <- dirname(dirname(revision_dir))
results_dir <- file.path(dirname(pilot_file), "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(revision_dir, "R", "manydist_bridge.R"))
source(file.path(revision_dir, "R", "distance_r1.R"))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("The pilot requires devtools to load the local manydist checkout.")
}
if (!requireNamespace("mclust", quietly = TRUE)) {
  stop("The pilot requires mclust for adjusted Rand index.")
}

devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)

pilot_reps <- as.integer(Sys.getenv("R1_PILOT_REPS", unset = "8"))
n_per_moon <- as.integer(Sys.getenv("R1_PILOT_N_PER_MOON", unset = "90"))
prop_nn_grid <- as.numeric(strsplit(
  Sys.getenv("R1_PILOT_PROP_NN", unset = "0.05,0.10,0.20"),
  ",",
  fixed = TRUE
)[[1L]])
if (!is.finite(pilot_reps) || pilot_reps < 2L) stop("R1_PILOT_REPS must be >= 2.")
if (!is.finite(n_per_moon) || n_per_moon < 18L || n_per_moon %% 9L != 0L) {
  stop("R1_PILOT_N_PER_MOON must be >= 18 and divisible by 9.")
}
if (!length(prop_nn_grid) || any(!is.finite(prop_nn_grid)) ||
    any(prop_nn_grid <= 0 | prop_nn_grid > 1)) {
  stop("R1_PILOT_PROP_NN must be a comma-separated grid in (0, 1].")
}
prop_nn_grid <- unique(prop_nn_grid)

make_moons <- function(seed, n_per_moon = 90L, interaction_signal = TRUE) {
  set.seed(seed)
  theta <- ((seq_len(n_per_moon) - 0.5) / n_per_moon) * pi
  noise_sd <- 0.055

  upper <- data.frame(
    V1 = cos(theta) + stats::rnorm(n_per_moon, sd = noise_sd),
    V2 = sin(theta) + stats::rnorm(n_per_moon, sd = noise_sd),
    arc = theta / pi,
    moon = "upper"
  )
  lower <- data.frame(
    V1 = 1 - cos(theta) + stats::rnorm(n_per_moon, sd = noise_sd),
    V2 = 0.50 - sin(theta) + stats::rnorm(n_per_moon, sd = noise_sd),
    arc = theta / pi,
    moon = "lower"
  )
  dat <- rbind(upper, lower)
  dat$moon <- factor(dat$moon, levels = c("upper", "lower"))
  equal_thirds <- rep(c("A", "B", "C"), each = n_per_moon / 3L)

  if (isTRUE(interaction_signal)) {
    # Each moon contains exactly one third A, B, and C along its arc.
    dat$band <- factor(
      rep(equal_thirds, times = 2L),
      levels = c("A", "B", "C")
    )
  } else {
    dat$band <- NA_character_
    for (lev in levels(dat$moon)) {
      idx <- which(dat$moon == lev)
      dat$band[idx] <- sample(equal_thirds)
    }
    dat$band <- factor(dat$band, levels = c("A", "B", "C"))
  }

  # Exactly balance nuisance within every moon-by-band cell, making nuisance
  # independent of both factors by construction, then randomize locations.
  dat$nuisance <- NA_character_
  for (moon_level in levels(dat$moon)) {
    for (band_level in levels(dat$band)) {
      idx <- which(dat$moon == moon_level & dat$band == band_level)
      stopifnot(length(idx) %% 3L == 0L)
      dat$nuisance[idx] <- sample(rep(
        c("N1", "N2", "N3"),
        each = length(idx) / 3L
      ))
    }
  }
  dat$nuisance <- factor(dat$nuisance, levels = c("N1", "N2", "N3"))

  truth <- if (isTRUE(interaction_signal)) {
    interaction(dat$moon, dat$band, sep = " x ", drop = TRUE)
  } else {
    dat$moon
  }

  list(
    x = dat[c("V1", "V2", "band", "nuisance")],
    truth = truth,
    scenario = if (interaction_signal) "interaction_moons" else "null_interaction"
  )
}

make_marginal_signal <- function(seed) {
  generated <- manydist::gen_mixed(
    k_true = 3,
    clustSizeEq = 60,
    numsignal = 3,
    numnoise = 2,
    catsignal = 3,
    catnoise = 2,
    q = 3,
    q_err = 4,
    numsep = 0.45,
    catsep = 0.45,
    seed = seed
  )
  list(
    x = generated$df[setdiff(names(generated$df), "y")],
    truth = generated$y,
    scenario = "marginal_signal_control"
  )
}

as_distance_matrix <- function(x) as.matrix(x$distance)

fit_distances <- function(dat, prop_nn, gate_seed) {
  r1_fit <- scmix_distance_r1(
    data = dat$x,
    prop_nn = prop_nn,
    score = "ba",
    decision = "prior_corrected",
    interaction_gate = "permutation_max",
    gate_permutations = 99L,
    gate_alpha = 0.01,
    gate_seed = gate_seed,
    strict_provenance = TRUE,
    manydist_repo = repo_root
  )

  legacy_int <- manydist::mdist(
    dat$x,
    preset = "custom",
    method_cat = "tvd",
    method_num = "pc_scores",
    commensurable = FALSE,
    interaction = TRUE,
    prop_nn = prop_nn,
    score = "ba",
    decision = "prior_corrected"
  )

  zero <- matrix(0, nrow(dat$x), nrow(dat$x))
  rho <- r1_fit$config$rho

  normalize_by_level_pairs <- function(observation_distance, level_delta) {
    denominator <- mean(level_delta[upper.tri(level_delta)])
    if (!is.finite(denominator) || denominator <= sqrt(.Machine$double.eps)) {
      return(zero)
    }
    observation_distance / denominator
  }

  # Diagnostic reconstruction of the submitted manuscript's separate
  # category-pair mean normalization. It is not the revised IFCS-talk method.
  submitted_normalized_components <- lapply(
    r1_fit$components$categorical,
    function(component) {
      tvd_normalized <- normalize_by_level_pairs(
        component$tvd_raw,
        component$tvd_level_delta
      )
      interaction_normalized <- normalize_by_level_pairs(
        component$interaction_raw,
        component$interaction_level_delta
      )
      (1 - rho) * tvd_normalized + rho * interaction_normalized
    }
  )
  submitted_normalized_diagnostic <- r1_fit$components$numeric$scaled +
    Reduce(`+`, submitted_normalized_components, init = zero)

  list(
    distances = list(
      AB_BW_R1_interaction = as.matrix(r1_fit$distance),
      AB_BW_R1_no_interaction = as.matrix(r1_fit$distance_no_interaction),
      submitted_normalized_diagnostic = submitted_normalized_diagnostic,
      u_dep_bw_package = as_distance_matrix(manydist::mdist(dat$x, preset = "u_dep_bw")),
      legacy_u_dep_interaction = as_distance_matrix(legacy_int),
      Gower = as_distance_matrix(manydist::mdist(dat$x, preset = "gower")),
      Euclidean_onehot = as_distance_matrix(manydist::mdist(dat$x, preset = "euclidean"))
    ),
    r1_fit = r1_fit
  )
}

cluster_and_score <- function(D, truth, affinity, seed) {
  set.seed(seed)
  predicted <- manydist::spectral_from_dist(
    D = D,
    k = nlevels(factor(truth)),
    affinity_method = affinity
  )
  mclust::adjustedRandIndex(predicted, truth)
}

result_rows <- list()
diagnostic_fits <- list()
row_id <- 0L

for (replicate_id in seq_len(pilot_reps)) {
  base_seed <- 20260918L + replicate_id * 100L
  datasets <- list(
    make_moons(base_seed + 1L, n_per_moon, interaction_signal = TRUE),
    make_moons(base_seed + 2L, n_per_moon, interaction_signal = FALSE),
    make_marginal_signal(base_seed + 3L)
  )

  for (dat in datasets) {
    for (prop_nn in prop_nn_grid) {
      fitted <- fit_distances(
        dat,
        prop_nn = prop_nn,
        gate_seed = base_seed +
          10000L +
          1000L * match(
            dat$scenario,
            c("interaction_moons", "null_interaction", "marginal_signal_control")
          ) +
          as.integer(round(100 * prop_nn))
      )
      diagnostic_key <- paste(
        dat$scenario,
        replicate_id,
        sprintf("nn%.2f", prop_nn),
        sep = "_"
      )
      diagnostic_fits[[diagnostic_key]] <- fitted$r1_fit

      for (affinity in c("gaussian", "selftune")) {
        for (method in names(fitted$distances)) {
          row_id <- row_id + 1L
          scored <- tryCatch(
            list(
              ARI = cluster_and_score(
                fitted$distances[[method]],
                dat$truth,
                affinity,
                # The same initialization stream is used for every method and
                # prop_nn value within a scenario/affinity replicate.
                seed = base_seed + match(affinity, c("gaussian", "selftune"))
              ),
              error = NA_character_
            ),
            error = function(e) {
              list(ARI = NA_real_, error = conditionMessage(e))
            }
          )
          result_rows[[row_id]] <- data.frame(
            scenario = dat$scenario,
            replicate = replicate_id,
            seed = base_seed,
            n = nrow(dat$x),
            prop_nn = prop_nn,
            method = method,
            affinity = affinity,
            ARI = scored$ARI,
            error = scored$error,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

results <- do.call(rbind, result_rows)

summarise_metric <- function(fun, output_name) {
  out <- aggregate(
    results["ARI"],
    by = results[c("scenario", "prop_nn", "affinity", "method")],
    FUN = function(x) fun(x[is.finite(x)])
  )
  names(out)[names(out) == "ARI"] <- output_name
  out
}

summary_results <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = c("scenario", "prop_nn", "affinity", "method"),
    all = TRUE
  ),
  list(
    summarise_metric(function(x) if (length(x)) mean(x) else NA_real_, "mean_ARI"),
    summarise_metric(
      function(x) if (length(x)) stats::median(x) else NA_real_,
      "median_ARI"
    ),
    summarise_metric(
      function(x) if (length(x) > 1L) stats::sd(x) else NA_real_,
      "sd_ARI"
    ),
    summarise_metric(length, "successful_replicates")
  )
)
summary_results$failed_replicates <- pilot_reps - summary_results$successful_replicates
summary_results <- summary_results[order(
  summary_results$scenario,
  summary_results$prop_nn,
  summary_results$affinity,
  -summary_results$mean_ARI
), ]

write.csv(results, file.path(results_dir, "pilot_results.csv"), row.names = FALSE)
write.csv(summary_results, file.path(results_dir, "pilot_summary.csv"), row.names = FALSE)
write.csv(
  results[!is.na(results$error), ],
  file.path(results_dir, "pilot_failures.csv"),
  row.names = FALSE
)
saveRDS(diagnostic_fits, file.path(results_dir, "pilot_diagnostics.rds"))
gate_diagnostics <- do.call(rbind, lapply(names(diagnostic_fits), function(key) {
  fit <- diagnostic_fits[[key]]
  data.frame(
    diagnostic_key = key,
    variable = names(fit$interaction_gate$detected),
    statistic = unname(fit$interaction_gate$statistic),
    p_value = unname(fit$interaction_gate$p_value),
    detected = unname(fit$interaction_gate$detected),
    permutations = fit$interaction_gate$permutations,
    alpha = fit$interaction_gate$alpha,
    seed = fit$interaction_gate$seed,
    stringsAsFactors = FALSE
  )
}))
write.csv(
  gate_diagnostics,
  file.path(results_dir, "pilot_gate_diagnostics.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))

cat("Pilot complete.\n")
cat("Replicates:", pilot_reps, "\n")
cat("Results:", file.path(results_dir, "pilot_results.csv"), "\n")
print(summary_results, row.names = FALSE)
