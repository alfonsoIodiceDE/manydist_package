#!/usr/bin/env Rscript

# Lightweight generator tuning: no interaction gate and no external
# competitors.  This checks whether one observation contains enough
# information to recover the categorical association regime before the
# confirmatory pilot is rerun.

pilot_dir <- normalizePath(
  "spectral_cluster_experiments/revision_r1/structured_blocks_pilot",
  mustWork = TRUE
)
revision_dir <- dirname(pilot_dir)
full_dir <- file.path(revision_dir, "full")
repo_root <- dirname(dirname(revision_dir))

devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)
source(file.path(full_dir, "R", "engine.R"))
source(file.path(pilot_dir, "R", "generator.R"))

base_specification <- yaml::read_yaml(file.path(pilot_dir, "design.yml"))$data
grid <- expand.grid(
  q_categorical = 24L,
  categorical_profile_strength = c(0.20, 0.35, 0.50),
  stringsAsFactors = FALSE
)
seeds <- c(26092501L, 26092601L)
rows <- list()
index <- 0L

cluster_distance <- function(distance, truth, seed) {
  affinity <- full_affinity(distance, "gaussian", list())
  fit <- full_spectral_cluster(
    affinity$matrix,
    k = nlevels(truth),
    seed = seed,
    kmeans_specification = list(
      nstart = 100L,
      iter_max = 100L,
      algorithm = "Hartigan-Wong"
    )
  )
  mclust::adjustedRandIndex(fit$cluster, truth)
}

for (grid_row in seq_len(nrow(grid))) {
  for (seed in seeds) {
    specification <- base_specification
    q <- grid$q_categorical[[grid_row]]
    specification$categorical_levels <- rep(c(2L, 3L, 4L), length.out = q)
    specification$categorical_model <- "nominal_association_blocks"
    specification$categorical_copy_probability <- 0.85
    specification$categorical_profile_strength <-
      grid$categorical_profile_strength[[grid_row]]
    generated <- structured_blocks_make_pair(specification, seed)

    for (arm in c("aligned", "decoupled")) {
      observed <- generated[[arm]]
      categorical <- observed$x[vapply(observed$x, is.factor, logical(1))]
      distance <- as.matrix(manydist::mdist(
        categorical,
        preset = "custom",
        method_cat = "tvd",
        method_num = "pc_scores",
        commensurable = FALSE,
        interaction = FALSE
      )$distance)
      index <- index + 1L
      rows[[index]] <- data.frame(
        q_categorical = q,
        categorical_profile_strength = specification$categorical_profile_strength,
        seed = seed,
        arm = arm,
        categorical_ARI = cluster_distance(
          distance,
          observed$truth,
          seed + 11L
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, rows)
print(results, digits = 4, row.names = FALSE)
print(
  stats::aggregate(
    categorical_ARI ~ q_categorical + categorical_profile_strength + arm,
    results,
    mean
  ),
  digits = 4,
  row.names = FALSE
)
