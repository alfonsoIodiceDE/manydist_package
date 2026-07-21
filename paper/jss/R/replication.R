run_replication <- function(root = paper_root()) {
  root <- normalizePath(root, mustWork = TRUE)
  ensure_output_dirs(root)
  package_version <- load_manydist_source(root)

  seeds <- paper_seeds()
  set.seed(unname(seeds[["distance"]]))

  penguins <- prepare_penguins()
  specifications <- paper_distance_specifications()
  validate_distance_specifications(specifications)
  distances <- compute_core_distances(
    penguins$predictors,
    specifications
  )

  write_csv_artifact(
    penguins_summary_table(penguins),
    paper_path("generated", "tables", "penguins-summary.csv", root = root)
  )
  write_csv_artifact(
    distance_specifications_table(specifications),
    paper_path(
      "generated",
      "tables",
      "distance-specifications.csv",
      root = root
    )
  )
  write_csv_artifact(
    core_distance_summary(distances),
    paper_path("generated", "tables", "core-distance-summary.csv", root = root)
  )
  write_csv_artifact(
    manydist::dist_methods_tbl(),
    paper_path("generated", "tables", "method-registry.csv", root = root)
  )
  write_csv_artifact(
    categorical_method_table(root),
    paper_path("generated", "tables", "categorical-methods.csv", root = root)
  )

  saveRDS(
    penguins,
    paper_path("generated", "results", "penguins.rds", root = root)
  )
  saveRDS(
    specifications,
    paper_path(
      "generated",
      "results",
      "distance-specifications.rds",
      root = root
    )
  )
  saveRDS(
    serializable_core_distances(distances),
    paper_path("generated", "results", "core-distances.rds", root = root),
    version = 3
  )

  write_mds_figure(
    distances$u_dep,
    penguins$response,
    paper_path(
      "generated",
      "figures",
      "fig-mds-unbiased-dependent.pdf",
      root = root
    )
  )

  session_lines <- c(
    paste0("manydist package version: ", package_version),
    paste0("repository root: ", repo_root(root)),
    "",
    capture.output(utils::sessionInfo())
  )
  writeLines(
    session_lines,
    paper_path("generated", "results", "session-info.txt", root = root)
  )

  manifest <- write_replication_manifest(root)
  message("Replication complete: ", manifest)
  manifest
}
