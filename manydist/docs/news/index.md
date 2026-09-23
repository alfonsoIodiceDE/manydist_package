# Changelog

## manydist 0.5.2

### Diagnostics

- Changed MDS-based diagnostics in
  [`lovo_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/lovo_mdist.md)
  and
  [`compare_lovo_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/compare_lovo_mdist.md)
  to be opt-in. Set `mds = TRUE` to compute congruence and alienation
  diagnostics; `dims = 2` remains the default MDS dimensionality when
  MDS is requested. Distance-based LOVO diagnostics remain available by
  default, and clustering diagnostics remain controlled by `cluster_k`.
- Added compact [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for
  `MDistBenchmark` objects. The full benchmark remains available as a
  tibble, while interactive output now emphasizes run status, failures,
  and pairwise diagnostic ranges.

### Distance construction

- Added the `"u_dep_bw"` preset for association-aware distances with
  block-wise numerical commensurability. The preset computes Manhattan
  distances on whitened principal-component scores and scales the
  complete numerical block so that its mean over distinct training pairs
  equals the number of original numerical variables. It retains the
  response-aware categorical construction used by `"u_dep"`.

- `"u_dep_bw"` supports principal-component selection through `ncomp` or
  `threshold`. For new observations, the PCA transformation and
  numerical block scaling estimated from the training data are reused.

## manydist 0.5.1

CRAN release: 2026-07-23

### Benchmarking and diagnostics

- Extended
  [`benchmark_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/benchmark_mdist.md)
  to compare every pair of successful distance specifications using mean
  absolute distance differences, symmetric relative distance,
  multidimensional-scaling congruence, and alienation.
- Added optional clustering comparisons to
  [`benchmark_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/benchmark_mdist.md).
  Supplying `cluster_k` computes pairwise adjusted Rand indices for PAM,
  hierarchical, and/or spectral clustering; clustering is skipped when
  `cluster_k = NULL`.
- Added
  [`benchmark_comparisons()`](https://alfonsoiodicede.github.io/manydist_package/reference/benchmark_comparisons.md)
  to extract the pairwise diagnostics stored in an `MDistBenchmark`
  result without recomputing the distances.
- Added an
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  method for `MDistBenchmark` objects, with annotated heatmaps for
  distance, geometry, and clustering-agreement diagnostics.

### Distance construction and recipe workflows

- Updated
  [`step_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/step_mdist.md)
  so response-aware specifications can obtain a single outcome directly
  from the recipe formula during preparation. The fitted response-aware
  profiles are reused when new data are baked, so assessment and test
  outcomes are neither required nor used.
- Added the `response_used` argument to
  [`step_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/step_mdist.md),
  allowing response use to be disabled explicitly.
- Allowed `method_num` to override the default standardization of the
  `"euclidean"` preset for numerical-only data. In particular,
  `method_num = "none"` computes ordinary Euclidean distances on the
  original variables.

### Data

- Added `wdi_2022`, a documented snapshot of selected 2022 World
  Development Indicators for reproducible mixed-type distance examples.

### Documentation and testing

- Added focused tests for response-aware
  [`step_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/step_mdist.md)
  workflows and the pairwise benchmarking interface.
- Expanded the package website with task-oriented articles on distance
  construction, diagnostics and benchmarking, clustering, and
  nearest-neighbour workflows.

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
  [`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md)
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
  [`step_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/step_mdist.md)
  for integrating `manydist` distances into `recipes` and tidymodels
  workflows.
- Added
  [`nearest_neighbor_dist()`](https://alfonsoiodicede.github.io/manydist_package/reference/nearest_neighbor_dist.md)
  and related prediction functions for nearest-neighbour models based on
  precomputed or manydist-generated distances.
- Added
  [`pam_dist()`](https://alfonsoiodicede.github.io/manydist_package/reference/pam_dist.md)
  for partitioning around medoids using manydist dissimilarities.
- Added
  [`spectral_dist()`](https://alfonsoiodicede.github.io/manydist_package/reference/spectral_dist.md)
  and
  [`spectral_from_dist()`](https://alfonsoiodicede.github.io/manydist_package/reference/spectral_from_dist.md)
  for spectral clustering from distance matrices.
- Added support functions for converting distances to affinities and
  fitting distance-based clustering models.

### Variable importance and diagnostics

- Added
  [`lovo_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/lovo_mdist.md)
  for leave-one-variable-out diagnostics of distance matrices.
- Added
  [`compare_lovo_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/compare_lovo_mdist.md)
  and
  [`lovo_method_spec()`](https://alfonsoiodicede.github.io/manydist_package/reference/lovo_method_spec.md)
  for comparing LOVO diagnostics across multiple distance
  specifications.
- Added congruence- and alienation-based diagnostics for comparing
  multidimensional scaling configurations.
- Added optional clustering-based LOVO diagnostics using PAM,
  hierarchical clustering, and spectral clustering.

### Data generation and benchmarking

- Added
  [`gen_mixed()`](https://alfonsoiodicede.github.io/manydist_package/reference/gen_mixed.md)
  and
  [`generate_dataset()`](https://alfonsoiodicede.github.io/manydist_package/reference/generate_dataset.md)
  for generating mixed-type example and simulation data.
- Added
  [`benchmark_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/benchmark_mdist.md)
  for benchmarking distance specifications across datasets and method
  grids.
- Added
  [`all_dist_method_specs()`](https://alfonsoiodicede.github.io/manydist_package/reference/all_dist_method_specs.md)
  and distance-method metadata helpers.

### Documentation

- Added documentation for the new modelling, clustering, LOVO,
  benchmarking, and recipe functions.
- Updated the package-level description and metadata for the `0.5.0`
  release.
- Removed CRAN-inappropriate development files, caches, and vignette
  outputs from the source build.
