# Higher-dimensional extension of the RSS26 interaction-moons example.
#
# The generator keeps the six-group latent truth fixed. Numerical variables
# are noisy smooth views of the same two-dimensional moon geometry. The
# categorical block contains one arc-band signal and up to four mutually
# pairwise-independent nuisance factors. A common row permutation of the
# complete categorical block creates the matched cross-type null without
# changing numerical data, truth, categorical margins, or categorical TVD.

.moon_scaling_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.moon_scaling_validate_integer <- function(x, name, lower = 1L) {
  if (length(x) != 1L || !is.finite(x) || x != as.integer(x) || x < lower) {
    stop("`", name, "` must be an integer >= ", lower, ".", call. = FALSE)
  }
  as.integer(x)
}

.moon_scaling_theta <- function(n_per_moon, probabilities) {
  counts <- as.integer(round(n_per_moon * probabilities))
  if (sum(counts) != n_per_moon || any(counts <= 0L)) {
    stop("Band probabilities do not give valid exact counts.", call. = FALSE)
  }
  boundaries <- c(0, cumsum(probabilities)) * pi
  theta <- unlist(lapply(seq_along(counts), function(j) {
    stats::runif(counts[[j]], boundaries[[j]], boundaries[[j + 1L]])
  }), use.names = FALSE)
  band <- factor(
    rep(c("A", "B", "C"), counts),
    levels = c("A", "B", "C")
  )
  list(theta = theta, band = band, counts = counts)
}

.moon_scaling_projection_angles <- function(p_max) {
  if (p_max == 2L) return(c(0, pi / 2))
  golden_fraction <- (sqrt(5) - 1) / 2
  c(0, pi / 2, ((seq_len(p_max - 2L) * golden_fraction) %% 1) * pi)
}

.moon_scaling_nuisance_bank <- function(moon, band, q_noise, seed) {
  if (q_noise == 0L) return(data.frame(row.names = seq_along(moon)))
  if (q_noise > 4L) {
    stop("The current exact ternary construction supports at most four nuisance factors.",
         call. = FALSE)
  }

  # OA(9, 4, 3, 2): every pair of columns contains every ordered pair once.
  oa <- expand.grid(u = 0:2, v = 0:2)
  oa <- cbind(
    oa$u,
    oa$v,
    (oa$u + oa$v) %% 3,
    (oa$u + 2 * oa$v) %% 3
  )
  levels_noise <- c("N1", "N2", "N3")
  truth_stratum <- interaction(moon, band, drop = TRUE)
  output <- matrix(NA_integer_, nrow = length(moon), ncol = q_noise)

  .moon_scaling_with_seed(seed, {
    for (stratum in levels(truth_stratum)) {
      rows <- which(truth_stratum == stratum)
      if (length(rows) %% nrow(oa) != 0L) {
        stop(
          "Each moon-by-band stratum must be divisible by 9 for exact nuisance balancing.",
          call. = FALSE
        )
      }
      block <- oa[rep(seq_len(nrow(oa)), length(rows) / nrow(oa)),
                  seq_len(q_noise), drop = FALSE]
      output[rows, ] <- block[sample.int(nrow(block)), , drop = FALSE]
    }
  })

  output <- as.data.frame(output, stringsAsFactors = TRUE)
  names(output) <- paste0("nuisance_", seq_len(q_noise))
  output[] <- lapply(output, function(x) factor(levels_noise[x + 1L], levels = levels_noise))
  output
}

moon_scaling_generate <- function(
    seed,
    n_per_moon = 180L,
    p_numeric = 2L,
    q_categorical = 2L,
    interaction_state = c("signal", "null"),
    moon_noise_sd = 0.08,
    sensor_noise_sd = 0.02,
    band_probabilities = c(A = 0.10, B = 0.45, C = 0.45),
    p_numeric_max = 10L,
    q_categorical_max = 5L
) {
  interaction_state <- match.arg(interaction_state)
  seed <- .moon_scaling_validate_integer(seed, "seed", lower = 0L)
  n_per_moon <- .moon_scaling_validate_integer(n_per_moon, "n_per_moon", lower = 18L)
  p_numeric <- .moon_scaling_validate_integer(p_numeric, "p_numeric", lower = 2L)
  q_categorical <- .moon_scaling_validate_integer(q_categorical, "q_categorical", lower = 2L)
  p_numeric_max <- .moon_scaling_validate_integer(p_numeric_max, "p_numeric_max", lower = 2L)
  q_categorical_max <- .moon_scaling_validate_integer(q_categorical_max, "q_categorical_max", lower = 2L)
  if (p_numeric > p_numeric_max || q_categorical > q_categorical_max) {
    stop("Requested dimensions exceed the declared maximum banks.", call. = FALSE)
  }
  if (q_categorical_max > 5L) {
    stop("`q_categorical_max` is currently limited to 5.", call. = FALSE)
  }
  if (length(moon_noise_sd) != 1L || !is.finite(moon_noise_sd) || moon_noise_sd < 0 ||
      length(sensor_noise_sd) != 1L || !is.finite(sensor_noise_sd) || sensor_noise_sd < 0) {
    stop("Noise standard deviations must be finite and nonnegative.", call. = FALSE)
  }
  if (length(band_probabilities) != 3L || any(!is.finite(band_probabilities)) ||
      any(band_probabilities <= 0) ||
      abs(sum(band_probabilities) - 1) > sqrt(.Machine$double.eps)) {
    stop("`band_probabilities` must contain three positive values summing to one.",
         call. = FALSE)
  }
  names(band_probabilities) <- c("A", "B", "C")

  upper_draw <- .moon_scaling_with_seed(seed + 1L, {
    .moon_scaling_theta(n_per_moon, band_probabilities)
  })
  lower_draw <- .moon_scaling_with_seed(seed + 2L, {
    .moon_scaling_theta(n_per_moon, band_probabilities)
  })
  theta <- c(upper_draw$theta, lower_draw$theta)
  latent_band <- factor(
    c(as.character(upper_draw$band), as.character(lower_draw$band)),
    levels = c("A", "B", "C")
  )
  moon <- factor(
    rep(c("upper", "lower"), each = n_per_moon),
    levels = c("upper", "lower")
  )

  clean_coordinates <- cbind(
    V1 = ifelse(moon == "upper", cos(theta), 1 - cos(theta)),
    V2 = ifelse(moon == "upper", sin(theta), 0.50 - sin(theta))
  )
  latent_coordinates <- clean_coordinates + .moon_scaling_with_seed(seed + 3L, {
    matrix(stats::rnorm(length(clean_coordinates), sd = moon_noise_sd),
           nrow = nrow(clean_coordinates), ncol = 2L)
  })

  angles <- .moon_scaling_projection_angles(p_numeric_max)
  loading_matrix <- cbind(cos(angles), sin(angles))
  numeric_bank <- latent_coordinates %*% t(loading_matrix)
  if (p_numeric_max > 2L && sensor_noise_sd > 0) {
    numeric_bank[, 3:p_numeric_max] <- numeric_bank[, 3:p_numeric_max, drop = FALSE] +
      .moon_scaling_with_seed(seed + 4L, {
        matrix(
          stats::rnorm(nrow(numeric_bank) * (p_numeric_max - 2L), sd = sensor_noise_sd),
          nrow = nrow(numeric_bank), ncol = p_numeric_max - 2L
        )
      })
  }
  colnames(numeric_bank) <- paste0("V", seq_len(p_numeric_max))

  nuisance_bank <- .moon_scaling_nuisance_bank(
    moon = moon,
    band = latent_band,
    q_noise = q_categorical_max - 1L,
    seed = seed + 5L
  )
  categorical_bank <- data.frame(
    band = latent_band,
    nuisance_bank,
    check.names = FALSE
  )

  # Randomize observation order once, before constructing the paired null.
  row_permutation <- .moon_scaling_with_seed(seed + 6L, {
    sample.int(length(moon))
  })
  numeric_bank <- numeric_bank[row_permutation, , drop = FALSE]
  clean_coordinates <- clean_coordinates[row_permutation, , drop = FALSE]
  theta <- theta[row_permutation]
  moon <- moon[row_permutation]
  latent_band <- latent_band[row_permutation]
  categorical_bank <- categorical_bank[row_permutation, , drop = FALSE]
  rownames(numeric_bank) <- NULL
  rownames(clean_coordinates) <- NULL
  rownames(categorical_bank) <- NULL

  categorical_permutation <- seq_len(nrow(categorical_bank))
  if (identical(interaction_state, "null")) {
    categorical_permutation <- .moon_scaling_with_seed(seed + 7L, {
      sample.int(nrow(categorical_bank))
    })
    categorical_bank <- categorical_bank[categorical_permutation, , drop = FALSE]
    rownames(categorical_bank) <- NULL
  }

  truth <- interaction(moon, latent_band, sep = " x ", drop = TRUE)
  selected_numeric <- as.data.frame(
    numeric_bank[, seq_len(p_numeric), drop = FALSE],
    check.names = FALSE
  )
  selected_categorical <- categorical_bank[, seq_len(q_categorical), drop = FALSE]
  data <- data.frame(selected_numeric, selected_categorical, check.names = FALSE)
  rownames(data) <- NULL

  list(
    data = data,
    truth = truth,
    latent = data.frame(
      moon = moon,
      band = latent_band,
      theta = theta,
      clean_V1 = clean_coordinates[, 1L],
      clean_V2 = clean_coordinates[, 2L]
    ),
    diagnostics = list(
      generator = "rss26_moons_high_dimension_v1",
      seed = seed,
      interaction_state = interaction_state,
      n = nrow(data),
      n_per_moon = n_per_moon,
      p_numeric = p_numeric,
      q_categorical = q_categorical,
      p_numeric_max = p_numeric_max,
      q_categorical_max = q_categorical_max,
      moon_noise_sd = moon_noise_sd,
      sensor_noise_sd = sensor_noise_sd,
      band_probabilities = band_probabilities,
      band_counts_per_moon = upper_draw$counts,
      observed_band_agreement = mean(data$band == latent_band),
      categorical_row_permutation = categorical_permutation,
      projection_angles = angles[seq_len(p_numeric)]
    )
  )
}
