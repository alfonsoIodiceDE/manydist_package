#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
pilot_dir <- dirname(dirname(script))
repo_root <- dirname(dirname(dirname(pilot_dir)))
source(file.path(pilot_dir, "R", "historical_feb2026_distance.R"))

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to verify against historical source.")
}
devtools::load_all(file.path(repo_root, "manydist"), quiet = TRUE)

old <- new.env(parent = asNamespace("manydist"))
for (filename in c("z_preproc.R", "cat_custom_delta.R", "cat_delta.R",
                   "cdist.R", "ndist.R", "mdist.R")) {
  source_lines <- suppressWarnings(system2(
    "git", c("-C", repo_root, "show", paste0("79563a7:manydist/R/", filename)),
    stdout = TRUE, stderr = TRUE
  ))
  if (!is.null(attr(source_lines, "status"))) {
    stop("Could not load historical ", filename, ": ",
         paste(source_lines, collapse = "\n"), call. = FALSE)
  }
  eval(parse(text = source_lines), envir = old)
}

set.seed(79563L)
numeric_data <- matrix(stats::rnorm(48L * 6L), nrow = 48L, ncol = 6L)
colnames(numeric_data) <- paste0("V", 1:6)
data <- data.frame(
  numeric_data,
  # Unequal category counts are essential here: a balanced test cannot
  # distinguish row-wise from column-wise normalization of the profile.
  C1 = factor(c(rep(1L, 8L), rep(2L, 15L), rep(3L, 25L))),
  C2 = factor(c(rep(1L, 31L), rep(2L, 17L))),
  C3 = factor(c(rep(1L, 11L), rep(2L, 37L))),
  C4 = factor(c(rep(3L, 7L), rep(2L, 18L), rep(1L, 23L)))
)
rebuilt <- historical_feb2026_sim2_distance(data, components = TRUE)
old_numeric <- old$ndist(
  x = data[1:6], method = "euclidean", commensurable = FALSE,
  scaling = "pc_scores", ncomp = ncol(data)
)
old_categorical <- old$cdist(
  x = data[7:10], method = "tot_var_dist", commensurable = FALSE
)$distance_mat
stopifnot(isTRUE(all.equal(rebuilt$numeric, as.matrix(old_numeric), tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(
  unname(rebuilt$categorical), unname(as.matrix(old_categorical)), tolerance = 1e-8
)))
stopifnot(isTRUE(all.equal(
  unname(rebuilt$distance),
  unname(as.matrix(old_numeric) + as.matrix(old_categorical)),
  tolerance = 1e-8
)))

# In the historical custom Euclidean branch, interaction never enters the
# calculation. Test the actual mdist() implementation, not only our formula.
old_without <- old$mdist(
  data, preset = "custom", distance_cont = "euclidean",
  distance_cat = "tot_var_dist", commensurable = FALSE,
  scaling_cont = "pc_scores", interaction = FALSE
)
old_with <- old$mdist(
  data, preset = "custom", distance_cont = "euclidean",
  distance_cat = "tot_var_dist", commensurable = FALSE,
  scaling_cont = "pc_scores", interaction = TRUE
)
stopifnot(identical(as.matrix(old_without$distance), as.matrix(old_with$distance)))
stopifnot(isTRUE(all.equal(
  unname(rebuilt$distance), unname(as.matrix(old_with$distance)), tolerance = 1e-8
)))
cat("Historical February comparator verification passed; interaction TRUE and FALSE are identical.\n")
