# Changelog

## manydist 0.5.0

CRAN release: 2026-06-09

### Major changes

- Expanded `manydist` from a package focused on mixed-type distance
  construction to a broader framework for distance-based learning with
  mixed-type data.
- Updated the package title and description to reflect support for
  distance construction, distance-based modelling workflows,
  variable-importance diagnostics, and clustering.
- Changed the package maintainer from Angelos Markos to Alfonso Iodice
  D’Enza.

### Distance construction

- Added a revised
  [`mdist()`](https://alfonsoiodicedenza.github.io/manydist/reference/mdist.md)
  interface and documentation for mixed-type distance construction.
- Added support for additional mixed-type distance specifications and
  presets.
- Added response-aware distance construction tools for supervised
  mixed-type workflows.
- Added interaction-aware distance components for continuous-categorical
  relationships.
- Added helper infrastructure for preprocessing and applying mixed-type
  distance specifications consistently across training and new data.
- Added utilities for generating and benchmarking distance-method
  specifications.

### Distance-based learning workflows

- Added
  [`step_mdist()`](https://alfonsoiodicedenza.github.io/manydist/reference/step_mdist.md)
  for integrating `manydist` distances into `recipes` and tidymodels
  workflows.
- Added
  [`nearest_neighbor_dist()`](https://alfonsoiodicedenza.github.io/manydist/reference/nearest_neighbor_dist.md)
  and related prediction functions for nearest-neighbour models based on
  precomputed or manydist-generated distances.
- Added
  [`pam_dist()`](https://alfonsoiodicedenza.github.io/manydist/reference/pam_dist.md)
  for partitioning around medoids using manydist dissimilarities.
- Added
  [`spectral_dist()`](https://alfonsoiodicedenza.github.io/manydist/reference/spectral_dist.md)
  and
  [`spectral_from_dist()`](https://alfonsoiodicedenza.github.io/manydist/reference/spectral_from_dist.md)
  for spectral clustering from distance matrices.
- Added support functions for converting distances to affinities and
  fitting distance-based clustering models.

### Variable importance and diagnostics

- Added
  [`lovo_mdist()`](https://alfonsoiodicedenza.github.io/manydist/reference/lovo_mdist.md)
  for leave-one-variable-out diagnostics of distance matrices.
- Added
  [`compare_lovo_mdist()`](https://alfonsoiodicedenza.github.io/manydist/reference/compare_lovo_mdist.md)
  and
  [`lovo_method_spec()`](https://alfonsoiodicedenza.github.io/manydist/reference/lovo_method_spec.md)
  for comparing LOVO diagnostics across multiple distance
  specifications.
- Added congruence- and alienation-based diagnostics for comparing
  multidimensional scaling configurations.
- Added optional clustering-based LOVO diagnostics using PAM,
  hierarchical clustering, and spectral clustering.

### Data generation and benchmarking

- Added
  [`gen_mixed()`](https://alfonsoiodicedenza.github.io/manydist/reference/gen_mixed.md)
  and
  [`generate_dataset()`](https://alfonsoiodicedenza.github.io/manydist/reference/generate_dataset.md)
  for generating mixed-type example and simulation data.
- Added
  [`benchmark_mdist()`](https://alfonsoiodicedenza.github.io/manydist/reference/benchmark_mdist.md)
  for benchmarking distance specifications across datasets and method
  grids.
- Added
  [`all_dist_method_specs()`](https://alfonsoiodicedenza.github.io/manydist/reference/all_dist_method_specs.md)
  and distance-method metadata helpers.

### Documentation

- Added documentation for the new modelling, clustering, LOVO,
  benchmarking, and recipe functions.
- Updated the package-level description and metadata for the `0.5.0`
  release.
- Removed CRAN-inappropriate development files, caches, and vignette
  outputs from the source build.
