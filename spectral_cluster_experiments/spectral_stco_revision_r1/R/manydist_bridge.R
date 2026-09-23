# Thin bridge between the paper-specific R1 analysis and manydist.
#
# Keep all use of non-exported manydist functions in this file.  The paper
# layer deliberately delegates the numerical u_dep_bw construction, the
# association-aware TVD category deltas, and the KNN interaction deltas to
# the package implementation pinned below.

SCMIX_EXPECTED_MANYDIST_VERSION <- "0.5.2"
SCMIX_EXPECTED_MANYDIST_COMMIT <- "ab27241f1a0846bb13087f804e23fd52746d1578"

.scmix_find_manydist_repo <- function(start = getwd()) {
  candidate <- normalizePath(start, mustWork = FALSE)

  repeat {
    if (file.exists(file.path(candidate, ".git")) &&
        file.exists(file.path(candidate, "manydist", "DESCRIPTION"))) {
      return(candidate)
    }

    parent <- dirname(candidate)
    if (identical(parent, candidate)) {
      return(NULL)
    }
    candidate <- parent
  }
}

.scmix_git_output <- function(repo, args) {
  if (is.null(repo) || !nzchar(repo) || !dir.exists(repo)) {
    return(character())
  }

  result <- suppressWarnings(
    system2(
      "git",
      c("-C", shQuote(normalizePath(repo)), args),
      stdout = TRUE,
      stderr = TRUE
    )
  )

  if (!is.null(attr(result, "status")) && attr(result, "status") != 0L) {
    return(character())
  }
  result
}

scmix_manydist_provenance <- function(
    manydist_repo = NULL,
    strict = TRUE,
    expected_version = SCMIX_EXPECTED_MANYDIST_VERSION,
    expected_commit = SCMIX_EXPECTED_MANYDIST_COMMIT
) {
  if (!requireNamespace("manydist", quietly = TRUE)) {
    stop("The paper-specific distance requires the `manydist` package.", call. = FALSE)
  }

  installed_version <- as.character(utils::packageVersion("manydist"))
  if (!identical(installed_version, expected_version)) {
    message <- paste0(
      "Expected manydist ", expected_version, " but found ", installed_version,
      ". Install or load the pinned package before running the revision."
    )
    if (isTRUE(strict)) stop(message, call. = FALSE) else warning(message, call. = FALSE)
  }

  repo <- manydist_repo
  if (is.null(repo)) repo <- .scmix_find_manydist_repo()
  if (!is.null(repo)) repo <- normalizePath(repo, mustWork = FALSE)

  head <- .scmix_git_output(repo, c("rev-parse", "HEAD"))
  head <- if (length(head)) head[[1L]] else NA_character_

  if (!is.na(head) && !identical(head, expected_commit)) {
    message <- paste0(
      "Expected manydist source commit ", expected_commit, " but found ", head,
      ". Use the pinned commit or set `strict_provenance = FALSE` for exploratory work."
    )
    if (isTRUE(strict)) stop(message, call. = FALSE) else warning(message, call. = FALSE)
  }

  source_status <- .scmix_git_output(
    repo,
    c("status", "--porcelain", "--", "manydist")
  )

  list(
    package_version = installed_version,
    package_path = system.file(package = "manydist"),
    expected_version = expected_version,
    source_repository = repo %||% NA_character_,
    source_commit = head,
    expected_commit = expected_commit,
    source_tree_dirty = length(source_status) > 0L,
    source_status = source_status,
    internal_functions = c("cat_delta", "delta_int_knn")
  )
}

.scmix_numeric_u_dep_bw <- function(x, ncomp = NULL, threshold = NULL) {
  args <- list(x = x, preset = "u_dep_bw")
  if (!is.null(ncomp)) args$ncomp <- ncomp
  if (!is.null(threshold)) args$threshold <- threshold

  fit <- do.call(manydist::mdist, args)
  scaled <- as.matrix(fit$distance)
  block_scale <- fit$params$numeric_block_scale

  if (is.null(block_scale) || !is.finite(block_scale) || block_scale <= 0) {
    stop("manydist did not return a valid u_dep_bw numerical block scale.", call. = FALSE)
  }

  list(
    fit = fit,
    raw = scaled / block_scale,
    scaled = scaled,
    block_mean = fit$params$numeric_block_mean,
    block_scale = block_scale,
    retained_ncomp = fit$params$retained_ncomp,
    preprocessor = fit$params$preprocessor$u_dep_bw_numeric
  )
}

.scmix_tvd_delta <- function(cat_data) {
  fn <- getFromNamespace("cat_delta", "manydist")
  fn(x = cat_data, method_cat = "tvd")
}

.scmix_interaction_delta <- function(
    numeric_distance,
    labels,
    prop_nn,
    score,
    decision,
    nn_order = NULL
) {
  fn <- getFromNamespace("delta_int_knn", "manydist")
  fn(
    D = numeric_distance,
    labels = labels,
    pi_nn = prop_nn,
    decision = decision,
    nn_order = nn_order,
    score = score
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
