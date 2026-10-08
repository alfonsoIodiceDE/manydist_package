# manydist 0.5.2

## Consistent new-data distances

- Replaced the Euclidean cross-distance backend used by `HLeucl` and `hl`
  with direct coordinate differences. This prevents zero test-to-training
  distances observed with larger training sets in the previous backend,
  for both commensurable and non-commensurable `HLeucl`.
- Corrected Gower's test-to-training average to divide by the number of
  original predictors rather than the number of expanded dummy columns.
  Training-to-training distances and the unaveraged Gower sum are unchanged.

## Training-only commensurability

- Corrected numerical and categorical commensurability so each contribution
  uses its training-to-training mean before aggregation, never the mean of a
  test-to-training matrix. This applies to `mdist()` and `step_mdist()`,
  including indicator-based categorical methods. Existing variable-wise
  training averages (including the zero diagonal) are preserved.
- Numerical absolute-distance means and categorical means are obtained without
  allocating observation-level training squares solely for normalization.
- Robust numerical preprocessing now also uses training medians and IQRs
  rather than recomputing them on each test batch.

## Direct prediction

- Added `knn_dist()` for direct classification, class probabilities, and
  regression from test-to-training distances, or from tabular predictors with
  a supplied distance function such as `mdist()`. The `response` argument
  selects and excludes the outcome column without manual predictor selection.
  Existing fit/predict engine
  functions remain available for tidymodels and compatibility.
- Fixed neighbour indexing for a single test observation with multiple
  neighbours and probability output for a single outcome level. Empty
  precomputed test matrices are also supported in the kNN helpers.

## Diagnostics

- Removed the redundant `spec_type` field from distance-specification grids
  and benchmark results. `preset = "custom"` selects explicit method settings;
  named presets select their predefined settings. Legacy `spec_type` columns
  are ignored on input. Update grid filters to use `preset` instead.

- Changed LOVO and LOVO-comparison plots so `reorder = TRUE` ranks variables
  within their categorical and numerical groups, preserving the meaning of
  the background bands. Comparison rankings use the mean metric across
  methods; `top_n` selection remains global.

- Changed the default fill gradient for pairwise benchmark heatmaps to
  coral (`#E76F51`) for lower values and blue (`#008CFF`) for higher values,
  with regular-weight white cell labels. This applies to all metrics plotted by
  `autoplot.MDistBenchmark()`, including clustering agreement.
- Left self-comparison tiles on the main diagonal blank. ARI heatmaps now
  use a fixed 0 to 1 colour scale instead of the observed range. Negative ARIs
  remain labelled and use the low-end colour.

- Changed MDS-based diagnostics in `lovo_mdist()` and
  `compare_lovo_mdist()` to be opt-in. Set `mds = TRUE` to compute
  congruence and alienation diagnostics; `dims = 2` remains the default MDS
  dimensionality when MDS is requested. Distance-based LOVO diagnostics remain
  available by default, and clustering diagnostics remain controlled by
  `cluster_k`.
- Added compact `print()` and `summary()` methods for `MDistBenchmark`
  objects. The full benchmark remains available as a tibble, while interactive
  output emphasizes run status and failures. `summary()` displays and invisibly
  returns the complete pairwise-results tibble.
- Removed `benchmark_comparisons()`. Use `pairs <- summary(result)` instead.
  This is a breaking API change; `autoplot()` continues to work directly on
  the benchmark object.

## Distance construction

- Added the `"u_dep_bw"` preset for association-aware distances with block-wise
  numerical commensurability. The preset computes Manhattan distances on
  whitened principal-component scores and scales the complete numerical block
  so that its mean over distinct training pairs equals the number of original
  numerical variables. It retains the response-aware categorical construction
  used by `"u_dep"`.

- `"u_dep_bw"` supports principal-component selection through `ncomp` or
  `threshold`. For new observations, the PCA transformation and numerical block
  scaling estimated from the training data are reused.

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
