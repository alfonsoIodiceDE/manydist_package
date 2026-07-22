# manydist 0.5.1

## Data

- Added `wdi_2022`, a documented snapshot of selected 2022 World Development
  Indicators for reproducible mixed-type distance examples.

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
