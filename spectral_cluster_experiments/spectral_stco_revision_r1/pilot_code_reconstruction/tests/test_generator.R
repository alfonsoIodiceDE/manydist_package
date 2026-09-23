#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
pilot_dir <- dirname(dirname(script))
source(file.path(pilot_dir, "R", "sim2_generator.R"))

for (rule in c("literal", "corrected")) {
  generated <- sim2_reconstructed_data(26091901L, n = 500L, index_rule = rule)
  repeated <- sim2_reconstructed_data(26091901L, n = 500L, index_rule = rule)
  stopifnot(identical(generated, repeated))
  stopifnot(identical(dim(generated$x), c(500L, 10L)))
  stopifnot(identical(as.integer(table(generated$x$C2)), c(400L, 100L)))
  stopifnot(identical(as.integer(table(generated$x$C3)), c(300L, 200L)))
  stopifnot(all(vapply(generated$x[7:10], is.factor, logical(1))))
  beta <- generated$diagnostics$coefficient_matrix[cbind(
    as.integer(generated$x$C1),
    as.integer(generated$x$C2)
  )]
  stopifnot(identical(
    as.character(generated$truth),
    ifelse(beta == 0.8, "C1", ifelse(beta == -0.4, "C2", "C3"))
  ))
  for (variable in 3:1) {
    expected <- ifelse(beta == 0.5, 8, 0) + beta *
      rowSums(generated$x[, paste0("V", (variable + 1L):6L), drop = FALSE])
    stopifnot(isTRUE(all.equal(generated$x[[paste0("V", variable)]], expected)))
  }
}

literal <- sim2_reconstructed_data(26091901L, 500L, "literal")
corrected <- sim2_reconstructed_data(26091901L, 500L, "corrected")
stopifnot(identical(as.integer(table(literal$x$C1)), c(200L, 59L, 241L)))
stopifnot(identical(as.integer(table(literal$x$C4)), c(200L, 89L, 211L)))
stopifnot(identical(as.integer(table(corrected$x$C1)), c(200L, 100L, 200L)))
stopifnot(identical(as.integer(table(corrected$x$C4)), c(200L, 150L, 150L)))

cat("Reconstructed Simulation 2 generator tests passed.\n")
