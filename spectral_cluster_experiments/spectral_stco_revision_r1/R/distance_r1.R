# Paper-specific AB-BW distance for the first revision.
#
# Source `manydist_bridge.R` before this file.  The reusable primitives stay
# in manydist; this file implements only the paper-specific nested composition,
# permutation gate, diagnostics, and provenance record.

.scmix_mean_distinct <- function(x) {
  stopifnot(is.matrix(x), nrow(x) == ncol(x))
  if (nrow(x) < 2L) return(NA_real_)
  mean(x[upper.tri(x)])
}

.scmix_validate_distance <- function(x, label) {
  if (!is.matrix(x) || nrow(x) != ncol(x)) {
    stop(label, " must be a square matrix.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(label, " contains non-finite values.", call. = FALSE)
  }
  if (any(x < -sqrt(.Machine$double.eps))) {
    stop(label, " contains negative values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-10))) {
    stop(label, " is not symmetric.", call. = FALSE)
  }

  x <- (x + t(x)) / 2
  x[x < 0] <- 0
  diag(x) <- 0
  x
}

.scmix_normalize_component <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  x <- .scmix_validate_distance(x, "Component dissimilarity")
  component_mean <- .scmix_mean_distinct(x)
  available <- is.finite(component_mean) && component_mean > tolerance

  if (available) {
    scale <- 1 / component_mean
    scaled <- x * scale
  } else {
    scale <- 0
    scaled <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  }

  list(
    raw = x,
    scaled = scaled,
    mean_distinct = component_mean,
    scale = scale,
    available = available
  )
}

.scmix_validate_unit_interval <- function(x, label, tolerance = 1e-10) {
  x <- .scmix_validate_distance(x, label)
  if (any(x > 1 + tolerance)) {
    stop(label, " must lie in [0, 1].", call. = FALSE)
  }
  x[x > 1] <- 1
  x
}

.scmix_column_names <- function(data, cols, predicate, label) {
  if (is.null(cols)) {
    return(names(data)[vapply(data, predicate, logical(1))])
  }

  if (is.numeric(cols)) {
    if (any(cols < 1L | cols > ncol(data)) || any(cols != as.integer(cols))) {
      stop(label, " contains invalid column positions.", call. = FALSE)
    }
    cols <- names(data)[as.integer(cols)]
  }

  if (!is.character(cols) || anyNA(cols) || any(!nzchar(cols))) {
    stop(label, " must contain column names or positions.", call. = FALSE)
  }
  if (anyDuplicated(cols)) stop(label, " contains duplicated columns.", call. = FALSE)
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    stop(label, " contains unknown columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  cols
}

.scmix_nn_order <- function(distance) {
  n <- nrow(distance)
  order_matrix <- matrix(NA_integer_, n, n - 1L)
  for (i in seq_len(n)) {
    ordered <- order(distance[, i], na.last = NA)
    order_matrix[i, ] <- ordered[ordered != i]
  }
  order_matrix
}

.scmix_pairwise_k <- function(labels, prop_nn) {
  labels <- droplevels(factor(labels))
  lev <- levels(labels)
  out <- matrix(NA_integer_, length(lev), length(lev), dimnames = list(lev, lev))
  if (length(lev) < 2L) return(out)

  for (a in seq_len(length(lev) - 1L)) {
    for (b in seq.int(a + 1L, length(lev))) {
      m <- sum(labels %in% lev[c(a, b)])
      k <- max(1L, min(m - 1L, as.integer(floor(prop_nn * (m - 1L)))))
      out[a, b] <- out[b, a] <- k
    }
  }
  out
}

.scmix_observation_distance <- function(level_delta, labels) {
  labels <- droplevels(factor(labels))
  indices <- as.integer(labels)
  result <- level_delta[indices, indices, drop = FALSE]
  dimnames(result) <- list(names(labels), names(labels))
  result
}

.scmix_component_status <- function(
    tvd_available,
    interaction_available,
    interaction_detected
) {
  if (!interaction_available) return("interaction_empirically_zero")
  if (!interaction_detected) return("interaction_not_detected")
  if (tvd_available) "interaction_detected_with_tvd" else "interaction_detected_tvd_zero"
}

.scmix_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.scmix_permutation_max_gate <- function(
    numeric_distance,
    categorical,
    observed_level_deltas,
    prop_nn,
    score,
    decision,
    nn_order,
    permutations,
    alpha,
    seed
) {
  observed_statistics <- vapply(observed_level_deltas, function(delta) {
    if (nrow(delta) < 2L) return(0)
    max(delta[upper.tri(delta)])
  }, numeric(1))

  null_maxima <- .scmix_with_seed(seed, vapply(seq_len(permutations), function(b) {
    # A common row permutation preserves the complete categorical block while
    # breaking its alignment with the numerical geometry.
    permutation <- sample.int(nrow(numeric_distance))
    max(vapply(seq_along(categorical), function(j) {
      delta <- .scmix_interaction_delta(
        numeric_distance = numeric_distance,
        labels = categorical[[j]][permutation],
        prop_nn = prop_nn,
        score = score,
        decision = decision,
        nn_order = nn_order
      )
      if (nrow(delta) < 2L) 0 else max(delta[upper.tri(delta)])
    }, numeric(1)))
  }, numeric(1)))

  p_values <- vapply(observed_statistics, function(statistic) {
    (1 + sum(null_maxima >= statistic)) / (permutations + 1)
  }, numeric(1))
  names(observed_statistics) <- names(categorical)
  names(p_values) <- names(categorical)

  list(
    method = "joint-label-permutation max statistic",
    detected = stats::setNames(p_values <= alpha, names(categorical)),
    p_value = p_values,
    statistic = observed_statistics,
    null_maxima = null_maxima,
    permutations = permutations,
    alpha = alpha,
    seed = seed
  )
}

#' Construct the paper-specific AB-BW interaction dissimilarity
#'
#' The numerical block is delegated to `manydist::mdist(...,
#' preset = "u_dep_bw")`.  Per-variable categorical TVD and KNN interaction
#' deltas are delegated to the pinned manydist internals. Their raw bounded
#' values enter a base-preserving, permutation-gated interaction extension.
#'
#' @param data Data frame containing the analysis variables only, or use the
#'   column selectors below to exclude identifiers/outcomes.
#' @param numeric_cols Names or positions of numerical variables.  Defaults to
#'   all numeric columns.
#' @param categorical_cols Names or positions of nominal variables.  Defaults
#'   to factor and character columns.
#' @param rho Added-interaction weight in `[0, 1]`. The primary paper
#'   specification defaults to `1 / q`, where `q` is the number of categorical
#'   variables; explicit alternatives are accepted for the pre-specified
#'   sensitivity analysis.
#' @param prop_nn Proportion used by the manydist KNN interaction estimator.
#' @param ncomp Optional number of whitened PCs.  The default retains all
#'   components with positive variance through `u_dep_bw`.
#' @param threshold Optional cumulative variance threshold; mutually exclusive
#'   with `ncomp`.
#' @param score,decision Passed to the manydist KNN interaction estimator.
#' @param interaction_gate Either `"permutation_max"` (the primary R1
#'   specification) or `"none"`. The permutation gate uses a common row
#'   permutation of the categorical block and a maximum statistic over every
#'   variable and category pair.
#' @param gate_permutations Number of permutations for the interaction gate.
#' @param gate_alpha Family-wise level for the maximum-statistic gate.
#' @param gate_seed Integer seed used locally by the permutation gate. The
#'   caller's random-number state is restored on exit.
#' @param strict_provenance Stop if the installed package version or detected
#'   source commit does not match the pinned revision configuration.
#' @param manydist_repo Optional path to the manydist source repository.
#'
#' @return A list of class `scmix_distance_r1` containing the final distance,
#'   every raw/scaled component, diagnostics, configuration, and provenance.
scmix_distance_r1 <- function(
    data,
    numeric_cols = NULL,
    categorical_cols = NULL,
    rho = NULL,
    prop_nn = 0.10,
    ncomp = NULL,
    threshold = NULL,
    score = "ba",
    decision = "prior_corrected",
    interaction_gate = c("permutation_max", "none"),
    gate_permutations = 99L,
    gate_alpha = 0.01,
    gate_seed = 1L,
    strict_provenance = TRUE,
    manydist_repo = NULL
) {
  if (!exists("scmix_manydist_provenance", mode = "function")) {
    stop("Source `R/manydist_bridge.R` before `R/distance_r1.R`.", call. = FALSE)
  }

  data <- as.data.frame(data)
  if (nrow(data) < 3L) stop("`data` must contain at least three rows.", call. = FALSE)
  if (is.null(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("`data` must have unique, non-empty column names.", call. = FALSE)
  }

  numeric_cols <- .scmix_column_names(data, numeric_cols, is.numeric, "`numeric_cols`")
  categorical_cols <- .scmix_column_names(
    data,
    categorical_cols,
    function(x) is.factor(x) || is.character(x),
    "`categorical_cols`"
  )

  overlap <- intersect(numeric_cols, categorical_cols)
  if (length(overlap)) {
    stop("Columns cannot be both numerical and categorical: ", paste(overlap, collapse = ", "), call. = FALSE)
  }
  if (!length(numeric_cols) || !length(categorical_cols)) {
    stop("The R1 distance requires at least one numerical and one categorical variable.", call. = FALSE)
  }
  if (!is.null(ncomp) && !is.null(threshold)) {
    stop("Specify only one of `ncomp` and `threshold`.", call. = FALSE)
  }
  if (length(prop_nn) != 1L || !is.finite(prop_nn) || prop_nn <= 0 || prop_nn > 1) {
    stop("`prop_nn` must be a finite number in (0, 1].", call. = FALSE)
  }
  score <- match.arg(score, c("ba", "logloss"))
  decision <- match.arg(decision, c("prior_corrected", "posterior"))
  interaction_gate <- match.arg(interaction_gate)
  if (length(gate_alpha) != 1L || !is.finite(gate_alpha) ||
      gate_alpha <= 0 || gate_alpha >= 1) {
    stop("`gate_alpha` must be a finite number in (0, 1).", call. = FALSE)
  }
  if (length(gate_permutations) != 1L || !is.finite(gate_permutations) ||
      gate_permutations < 1L || gate_permutations != as.integer(gate_permutations)) {
    stop("`gate_permutations` must be a positive integer.", call. = FALSE)
  }
  gate_permutations <- as.integer(gate_permutations)
  if (identical(interaction_gate, "permutation_max") &&
      1 / (gate_permutations + 1) > gate_alpha) {
    stop(
      "`gate_permutations` is too small to attain `gate_alpha`; use at least ",
      ceiling(1 / gate_alpha) - 1L, ".",
      call. = FALSE
    )
  }
  if (length(gate_seed) != 1L || !is.finite(gate_seed) ||
      gate_seed < 0 || gate_seed != as.integer(gate_seed)) {
    stop("`gate_seed` must be a nonnegative integer.", call. = FALSE)
  }
  gate_seed <- as.integer(gate_seed)

  numerical <- data[numeric_cols]
  categorical <- data[categorical_cols]

  if (anyNA(numerical) || any(!vapply(numerical, function(x) all(is.finite(x)), logical(1)))) {
    stop("Numerical variables must be finite and complete for the pilot.", call. = FALSE)
  }
  if (any(vapply(categorical, is.ordered, logical(1)))) {
    stop("Ordered factors are outside the nominal-variable scope of this R1 implementation.", call. = FALSE)
  }
  categorical <- lapply(categorical, function(x) droplevels(factor(x)))
  categorical <- as.data.frame(categorical, stringsAsFactors = TRUE)
  if (anyNA(categorical)) stop("Categorical variables must be complete for the pilot.", call. = FALSE)
  low_levels <- names(categorical)[vapply(categorical, nlevels, integer(1)) < 2L]
  if (length(low_levels)) {
    stop("Categorical variables need at least two observed levels: ", paste(low_levels, collapse = ", "), call. = FALSE)
  }

  q <- ncol(categorical)
  default_rho <- 1 / q
  if (is.null(rho)) rho <- default_rho
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho > 1) {
    stop("`rho` must be a finite number in [0, 1].", call. = FALSE)
  }

  provenance <- scmix_manydist_provenance(
    manydist_repo = manydist_repo,
    strict = strict_provenance
  )

  numeric_component <- .scmix_numeric_u_dep_bw(
    numerical,
    ncomp = ncomp,
    threshold = threshold
  )
  numeric_component$raw <- .scmix_validate_distance(
    numeric_component$raw,
    "Raw numerical dissimilarity"
  )
  numeric_component$scaled <- .scmix_validate_distance(
    numeric_component$scaled,
    "Scaled numerical dissimilarity"
  )
  numeric_component$mean_distinct_raw <- .scmix_mean_distinct(numeric_component$raw)
  numeric_component$mean_distinct_scaled <- .scmix_mean_distinct(numeric_component$scaled)
  numeric_component$fit <- NULL

  n <- nrow(data)
  row_ids <- rownames(data)
  dimnames(numeric_component$raw) <- list(row_ids, row_ids)
  dimnames(numeric_component$scaled) <- list(row_ids, row_ids)

  # The TVD delta is jointly estimated across categorical variables. With a
  # single categorical variable there is no other categorical profile, so its
  # TVD component is zero.
  tvd_output <- if (q > 1L) .scmix_tvd_delta(categorical) else NULL
  tvd_full <- if (is.null(tvd_output)) NULL else as.matrix(tvd_output$tvd)
  level_counts <- vapply(categorical, nlevels, integer(1))
  level_starts <- c(1L, head(cumsum(level_counts), -1L) + 1L)
  level_stops <- cumsum(level_counts)

  nn_order <- .scmix_nn_order(numeric_component$scaled)
  categorical_components <- vector("list", q)
  names(categorical_components) <- names(categorical)

  for (j in seq_len(q)) {
    labels <- categorical[[j]]
    lev <- levels(labels)

    if (is.null(tvd_full)) {
      tvd_level <- matrix(0, length(lev), length(lev), dimnames = list(lev, lev))
    } else {
      idx <- level_starts[j]:level_stops[j]
      tvd_level <- tvd_full[idx, idx, drop = FALSE]
      dimnames(tvd_level) <- list(lev, lev)
    }
    tvd_level <- .scmix_validate_unit_interval(
      tvd_level,
      paste0("TVD level dissimilarity for `", names(categorical)[j], "`")
    )
    tvd_raw <- .scmix_observation_distance(tvd_level, labels)
    dimnames(tvd_raw) <- list(row_ids, row_ids)
    tvd <- .scmix_normalize_component(tvd_raw)

    interaction_level <- .scmix_interaction_delta(
      numeric_distance = numeric_component$scaled,
      labels = labels,
      prop_nn = prop_nn,
      score = score,
      decision = decision,
      nn_order = nn_order
    )
    interaction_level <- .scmix_validate_unit_interval(
      interaction_level,
      paste0("Interaction level dissimilarity for `", names(categorical)[j], "`")
    )
    interaction_raw <- .scmix_observation_distance(interaction_level, labels)
    dimnames(interaction_raw) <- list(row_ids, row_ids)
    interaction <- .scmix_normalize_component(interaction_raw)

    categorical_components[[j]] <- list(
      variable = names(categorical)[j],
      levels = lev,
      tvd_level_delta = tvd_level,
      interaction_level_delta = interaction_level,
      tvd_raw = tvd$raw,
      tvd_scaled = tvd$scaled,
      interaction_raw = interaction$raw,
      interaction_scaled = interaction$scaled,
      tvd_available = tvd$available,
      interaction_available = interaction$available,
      tvd_mean_distinct = tvd$mean_distinct,
      interaction_mean_distinct = interaction$mean_distinct,
      tvd_scale = tvd$scale,
      interaction_scale = interaction$scale,
      pairwise_k = .scmix_pairwise_k(labels, prop_nn)
    )
  }

  gate <- if (identical(interaction_gate, "permutation_max")) {
    .scmix_permutation_max_gate(
      numeric_distance = numeric_component$scaled,
      categorical = categorical,
      observed_level_deltas = lapply(
        categorical_components,
        `[[`,
        "interaction_level_delta"
      ),
      prop_nn = prop_nn,
      score = score,
      decision = decision,
      nn_order = nn_order,
      permutations = gate_permutations,
      alpha = gate_alpha,
      seed = gate_seed
    )
  } else {
    list(
      method = "none",
      detected = stats::setNames(rep(TRUE, q), names(categorical)),
      p_value = stats::setNames(rep(NA_real_, q), names(categorical)),
      statistic = stats::setNames(vapply(
        categorical_components,
        function(component) {
          max(component$interaction_level_delta[
            upper.tri(component$interaction_level_delta)
          ])
        },
        numeric(1)
      ), names(categorical)),
      null_maxima = numeric(),
      permutations = 0L,
      alpha = gate_alpha,
      seed = gate_seed
    )
  }

  for (j in seq_len(q)) {
    component <- categorical_components[[j]]
    detected <- isTRUE(gate$detected[[j]])
    gated_interaction <- if (detected) {
      component$interaction_raw
    } else {
      matrix(0, n, n, dimnames = list(row_ids, row_ids))
    }
    combined <- component$tvd_raw + rho * gated_interaction
    combined <- .scmix_validate_distance(
      combined,
      paste0("Combined categorical dissimilarity for `", names(categorical)[j], "`")
    )
    dimnames(combined) <- list(row_ids, row_ids)

    component$interaction_detected <- detected
    component$interaction_gate_p_value <- gate$p_value[[j]]
    component$interaction_gate_statistic <- gate$statistic[[j]]
    component$interaction_gated_raw <- gated_interaction
    component$combined <- combined
    component$requested_weights <- c(tvd = 1, interaction = rho)
    component$effective_weights <- c(
      tvd = 1,
      interaction = if (detected) rho else 0
    )
    component$availability_reason <- .scmix_component_status(
      component$tvd_available,
      component$interaction_available,
      detected
    )
    categorical_components[[j]] <- component
  }

  categorical_sum <- Reduce(
    `+`,
    lapply(categorical_components, `[[`, "combined"),
    init = matrix(0, n, n, dimnames = list(row_ids, row_ids))
  )
  categorical_tvd_sum <- Reduce(
    `+`,
    lapply(categorical_components, `[[`, "tvd_raw"),
    init = matrix(0, n, n, dimnames = list(row_ids, row_ids))
  )
  categorical_tvd_scaled_sum <- Reduce(
    `+`,
    lapply(categorical_components, `[[`, "tvd_scaled"),
    init = matrix(0, n, n, dimnames = list(row_ids, row_ids))
  )
  categorical_interaction_gated_sum <- Reduce(
    `+`,
    lapply(categorical_components, `[[`, "interaction_gated_raw"),
    init = matrix(0, n, n, dimnames = list(row_ids, row_ids))
  )
  final <- numeric_component$scaled + categorical_sum
  final <- .scmix_validate_distance(final, "Final R1 dissimilarity")
  dimnames(final) <- list(row_ids, row_ids)
  no_interaction <- numeric_component$scaled + categorical_tvd_sum
  no_interaction <- .scmix_validate_distance(
    no_interaction,
    "R1 dissimilarity without interaction"
  )
  dimnames(no_interaction) <- list(row_ids, row_ids)

  unused_cols <- setdiff(names(data), c(numeric_cols, categorical_cols))
  diagnostics <- list(
    n = n,
    p_numeric = length(numeric_cols),
    q_categorical = q,
    retained_ncomp = numeric_component$retained_ncomp,
    numeric_mean_distinct = numeric_component$mean_distinct_scaled,
    categorical_mean_distinct = vapply(
      categorical_components,
      function(x) .scmix_mean_distinct(x$combined),
      numeric(1)
    ),
    tvd_available = vapply(categorical_components, `[[`, logical(1), "tvd_available"),
    interaction_available = vapply(categorical_components, `[[`, logical(1), "interaction_available"),
    interaction_detected = vapply(categorical_components, `[[`, logical(1), "interaction_detected"),
    interaction_gate_p_value = vapply(
      categorical_components,
      `[[`,
      numeric(1),
      "interaction_gate_p_value"
    ),
    availability_reason = vapply(categorical_components, `[[`, character(1), "availability_reason"),
    final_mean_distinct = .scmix_mean_distinct(final),
    unused_columns = unused_cols
  )

  structure(
    list(
      distance = final,
      distance_no_interaction = no_interaction,
      components = list(
        numeric = numeric_component,
        categorical = categorical_components,
        categorical_sum = categorical_sum,
        categorical_tvd_sum = categorical_tvd_sum,
        categorical_tvd_scaled_sum = categorical_tvd_scaled_sum,
        categorical_interaction_gated_sum = categorical_interaction_gated_sum
      ),
      interaction_gate = gate,
      diagnostics = diagnostics,
      config = list(
        numeric_cols = numeric_cols,
        categorical_cols = categorical_cols,
        rho = rho,
        rho_default = default_rho,
        rho_is_default = isTRUE(all.equal(rho, default_rho)),
        prop_nn = prop_nn,
        ncomp = ncomp,
        threshold = threshold,
        score = score,
        decision = decision,
        interaction_gate = interaction_gate,
        gate_permutations = if (interaction_gate == "permutation_max") {
          gate_permutations
        } else {
          0L
        },
        gate_alpha = gate_alpha,
        gate_seed = gate_seed,
        categorical_scale_for_final = "raw_unit_interval",
        diagnostic_normalization = "mean_distinct_pairs_not_used_in_final",
        requested_blend = "TVD_raw + rho * gate * INT_raw",
        availability_rule = paste(
          "the no-interaction base is retained; only permutation-detected",
          "interaction components are added"
        )
      ),
      provenance = provenance
    ),
    class = "scmix_distance_r1"
  )
}

print.scmix_distance_r1 <- function(x, ...) {
  cat("Paper-specific AB-BW R1 dissimilarity\n")
  cat("  observations: ", x$diagnostics$n, "\n", sep = "")
  cat("  numerical variables: ", x$diagnostics$p_numeric, "\n", sep = "")
  cat("  categorical variables: ", x$diagnostics$q_categorical, "\n", sep = "")
  cat("  retained PCs: ", x$diagnostics$retained_ncomp, "\n", sep = "")
  cat("  rho: ", format(x$config$rho), "\n", sep = "")
  cat(
    "  detected interactions: ",
    sum(x$diagnostics$interaction_detected),
    "/",
    x$diagnostics$q_categorical,
    "\n",
    sep = ""
  )
  invisible(x)
}
