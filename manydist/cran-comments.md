## Test environments

- local macOS Tahoe 26.5.2, R 4.6.0
- win-builder, R-devel (2026-07-22 r90289)
- win-builder, R-release (R 4.6.1)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Update

This is an update to version 0.5.1.

Main changes:

- Extended `benchmark_mdist()` with pairwise distance, geometry, and optional
  clustering-agreement diagnostics.
- Added `benchmark_comparisons()` and an `autoplot()` method for
  `MDistBenchmark` objects.
- Updated `step_mdist()` so response-aware specifications can obtain the
  outcome from a recipe and reuse the fitted profiles when baking new data.
- Added explicit control over response use in `step_mdist()` and over numerical
  preprocessing for the Euclidean preset on numerical-only data.
- Added the documented `wdi_2022` dataset.
- Added tests and expanded documentation for the new interfaces.

The Quarto articles used by the package website are not included in the CRAN
source package.
