compute_core_distances <- function(predictors, specifications) {
  validate_distance_specifications(specifications)

  purrr::map(specifications, function(specification) {
    do.call(
      manydist::mdist,
      c(list(x = predictors), specification)
    )
  })
}

core_distance_summary <- function(distances) {
  purrr::imap_dfr(
    distances,
    function(object, name) {
      matrix <- as.matrix(object$distance)
      values <- matrix[lower.tri(matrix)]

      tibble::tibble(
        specification = name,
        preset = object$preset,
        observations = nrow(matrix),
        mean_distance = mean(values),
        sd_distance = stats::sd(values),
        min_distance = min(values),
        max_distance = max(values)
      )
    }
  )
}

serializable_core_distances <- function(distances) {
  purrr::map(distances, function(object) {
    stable_params <- object$params[
      setdiff(names(object$params), "preprocessor")
    ]

    list(
      distance = as.matrix(object$distance),
      preset = object$preset,
      params = stable_params
    )
  })
}

write_mds_figure <- function(distance, species, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  coordinates <- stats::cmdscale(distance$distance, k = 2)
  palette <- c("#D73027", "#1A9850", "#4575B4")
  species <- droplevels(species)
  colours <- palette[as.integer(species)]

  grDevices::pdf(
    path,
    width = 7.2,
    height = 5.2,
    useDingbats = FALSE,
    title = "MDS - unbiased dependent distance",
    author = "manydist paper replication",
    timestamp = FALSE,
    producer = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::plot(
    coordinates[, 1],
    coordinates[, 2],
    col = colours,
    pch = 19,
    cex = 0.7,
    xlab = "Dimension 1",
    ylab = "Dimension 2",
    main = "MDS - unbiased dependent distance"
  )
  graphics::legend(
    "topright",
    legend = levels(species),
    col = palette[seq_along(levels(species))],
    pch = 19,
    bty = "n"
  )

  invisible(normalizePath(path, mustWork = FALSE))
}
