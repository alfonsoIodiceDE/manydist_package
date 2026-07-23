# manydist 0.5.1

## Benchmarking and diagnostics

- Extended `benchmark_mdist()` to compare every pair of successful distance
  specifications using mean absolute distance differences, symmetric relative
  distance, multidimensional-scaling congruence, and alienation.
- Added optional clustering comparisons to `benchmark_mdist()`. Supplying
  `cluster_k` computes pairwise adjusted Rand indices for PAM, hierarchical,
  and/or spectral clustering; clustering is skipped when `cluster_k = NULL`.
- Added `benchmark_comparisons()` to extract the pairwise diagnostics stored in
  an `MDistBenchmark` result without recomputing the distances.
- Added an `autoplot()` method for `MDistBenchmark` objects, with annotated
  heatmaps for distance, geometry, and clustering-agreement diagnostics.

## Distance construction and recipe workflows

- Updated `step_mdist()` so response-aware specifications can obtain a single
  outcome directly from the recipe formula during preparation. The fitted
  response-aware profiles are reused when new data are baked, so assessment
  and test outcomes are neither required nor used.
- Added the `response_used` argument to `step_mdist()`, allowing response use to
  be disabled explicitly.
- Allowed `method_num` to override the default standardization of the
  `"euclidean"` preset for numerical-only data. In particular,
  `method_num = "none"` computes ordinary Euclidean distances on the original
  variables.

## Data

- Added `wdi_2022`, a documented snapshot of selected 2022 World Development
  Indicators for reproducible mixed-type distance examples.

## Documentation and testing

- Added focused tests for response-aware `step_mdist()` workflows and the
  pairwise benchmarking interface.
- Expanded the package website with task-oriented articles on distance
  construction, diagnostics and benchmarking, clustering, and
  nearest-neighbour workflows.

# manydist 0.5.0

## Major changes

- Expanded `manydist` from a package focused on mixed-type distance construction to a broader framework for distance-based learning with mixed-type data.
- Updated the package title and description to reflect support for distance construction, distance-based modelling workflows, variable-importance diagnostics, and clustering.
- Changed the package maintainer from Angelos Markos to Alfonso Iodice D'Enza.

## Distance construction

- Added a revised `mdist()` interface and documentation for mixed-type distance construction.
- Added support for additional mixed-type distance specifications and presets.
- Added response-aware distance construction tools for supervised mixed-type workflows.
- Added interaction-aware distance components for continuous-categorical relationships.
- Added helper infrastructure for preprocessing and applying mixed-type distance specifications consistently across training and new data.
- Added utilities for generating and benchmarking distance-method specifications.

## Distance-based learning workflows

- Added `step_mdist()` for integrating `manydist` distances into `recipes` and tidymodels workflows.
- Added `nearest_neighbor_dist()` and related prediction functions for nearest-neighbour models based on precomputed or manydist-generated distances.
- Added `pam_dist()` for partitioning around medoids using manydist dissimilarities.
- Added `spectral_dist()` and `spectral_from_dist()` for spectral clustering from distance matrices.
- Added support functions for converting distances to affinities and fitting distance-based clustering models.

## Variable importance and diagnostics

- Added `lovo_mdist()` for leave-one-variable-out diagnostics of distance matrices.
- Added `compare_lovo_mdist()` and `lovo_method_spec()` for comparing LOVO diagnostics across multiple distance specifications.
- Added congruence- and alienation-based diagnostics for comparing multidimensional scaling configurations.
- Added optional clustering-based LOVO diagnostics using PAM, hierarchical clustering, and spectral clustering.

## Data generation and benchmarking

- Added `gen_mixed()` and `generate_dataset()` for generating mixed-type example and simulation data.
- Added `benchmark_mdist()` for benchmarking distance specifications across datasets and method grids.
- Added `all_dist_method_specs()` and distance-method metadata helpers.

## Documentation

- Added documentation for the new modelling, clustering, LOVO, benchmarking, and recipe functions.
- Updated the package-level description and metadata for the `0.5.0` release.
- Removed CRAN-inappropriate development files, caches, and vignette outputs from the source build.
