# Review of the current JSS manuscript — 2026-10-08

## Overall assessment

The shorter Sections 3–5 are easier to follow and put the package interface at
the centre of the article. The transition from distance construction to direct
diagnostics, downstream diagnostics, and learning workflows is clearer. Removing
the extended PCA development is reasonable for a software paper, provided the
remaining description accurately identifies the implemented transformation.

This review checks the current manuscript against the local package source; it
is not a fresh audit of the cited literature. Substantive wording and distance
algorithms were left unchanged. Manuscript edits were limited to the
benchmark-specification paragraph needed for the agreed API simplification
and two colour-markup repairs needed to render equation references correctly.

## Implementation issues to resolve before submission

### 1. Commensurability weights are not always held fixed for test data

**Resolved on 2026-10-08 in the subsequent normalization fix.** Numerical,
categorical, and indicator-based contributions now use training-only means
before aggregation. Tests cover changes in test-batch composition, row order,
single-observation application, and `step_mdist()`. Robust numerical scaling
also now uses training medians and IQRs. The reproduction below describes the
former behavior, not the corrected implementation. The existing full-square
averaging convention is preserved, so issue 3 remains open. Issue 2 was
resolved in the subsequent new-data calculation fix below.

The manuscript says that training-estimated weights are reused unchanged for
test observations. However, `ndist()` divides each test-to-training contribution
by the mean of that rectangular matrix (`manydist/R/ndist.R`, lines 99–109).
The categorical helper likewise uses means of the rectangular matrices
(`manydist/R/commensurability_for_cat.R`, lines 58–63).

Reproduction: training values `c(0, 2, 8, 10)`, test values `c(1, 9)`, and
`preset = "u_indep"`. The first test-to-training distance is 0.2222222.
Appending a further test value of 100 changes that same distance to 0.02884615.
The training data and original test observation have not changed.

This verifies test-batch-dependent normalization, not a change in Gower's
training ranges. In a single-variable example the change is a common scale
factor; with several variable-specific weights it can also change relative
contributions. Fitting inside each resample does not by itself fix this issue,
because the shared application code still uses the evaluation batch's means.
Either the package must store and reuse training-only weights, or the paper's
stronger fit-and-apply claims must be narrowed. This is separate from the
requested specification-table change and has not been patched here.

### 2. Mixed-data Gower normalization differs between the two application paths

**Resolved on 2026-10-08.** The new-data path now divides by the number of
original predictors. Tests compare numerical-only, categorical-only, and mixed
data against `cluster::daisy()`, including applying the distance to its own
training data, individual test observations, and recipe output. The unaveraged
sum is unchanged. The reproduction below documents the former denominator.

Gower's numerical ranges do come only from training data. Nevertheless, the
test-to-training branch divides by the number of baked columns, including dummy
columns, rather than the number of original variables (`manydist/R/mdist.R`,
lines 666–681).

With training data `z = c(0, 2, 8, 10)` and a factor `g = c("a", "a", "b", "b")`,
ordinary training-to-training Gower gives first-row distances
`c(0, 0.1, 0.9, 1)`. Passing the identical data through `new_data` gives
`c(0, 0.04, 0.36, 0.4)`. These should agree for complete data with the same
specification. Here the discrepancy is a common scale factor, so it does not
change nearest-neighbour rankings; it still matters for consistent distance
values and the claim of a conventional Gower scale.

### 3. The empirical-mean convention does not match Equation @eq-empirical-weights

The manuscript averages over distinct pairs with denominator `I(I-1)`.
The ordinary numerical and categorical commensurability implementations use
the whole square matrix, including zero diagonal entries. For four observations,
the resulting mean contribution over distinct pairs is 4/3, not one.
The block-scaled numerical implementation instead uses distinct pairs through
`stats::dist()`. This requires an explicit, consistent convention, especially
when comparing numerical-block and categorical-block totals. The displayed
formula is mathematically sound; the problem is correspondence with the code.

### 4. HLeucl / hl new-data Euclidean backend — resolved

The broader non-commensurable audit reproduced zero cross-distances from the
previous `Rfast::dista()` backend on larger examples in this installation
(Rfast 2.1.5.2), while the training path returned nonzero Euclidean distances.
This affected both non-commensurable and commensurable `HLeucl`, as well as
the categorical part of the `hl` preset. The earlier four-row tests did not
trigger the backend failure.

Replaced both calls with a direct Euclidean cross-distance calculation based
on coordinate differences, retaining the same training-derived dummy
representation and weights. Removed the now-unused Rfast package dependency.
Tests use 36 training rows, independent Euclidean references, self-application,
single test observations, enlarged test batches, and recipes. Training-distance
construction is unchanged.

After these fixes, 71 of the 72 specifications in the broader audit pass all
checks (batch independence, individual predictions, row order, self-application,
and recipe consistency). `kulczynski_s` rejects the audit's sparse categorical
profiles with a non-finite-value error; it is not a transformation-leakage
failure. `dkss`, `gudmm`, and `mod_gower` explicitly reject new-data calls.
Independent numerical references also agree for scaling and PCA, including a
separate reduced-PCA check. On the paper's penguin split, the corrected `hl`
changes three of 84 class predictions (accuracy 81/84 to 84/84); Gower's
correction rescales distances without changing predictions.

## Manuscript points to revise or discuss

1. **Shortened association-aware paragraph (Section 3).** The groupwise option
   currently implemented as `u_dep_bw` whitens the PC scores before computing
   Manhattan distances. Thus “preserves the relative importance of the
   principal components” should not suggest preservation of the original
   explained-variance weighting. A compact, accurate description would say it
   scales the complete Manhattan distance on whitened PC scores by one common
   factor. It need not restore the removed lengthy derivation.
2. **Complete preset table (Section 4).** `u_dep_bw` is available in the package
   but absent from the table described as complete. Add it or explicitly limit
   the table to the presets discussed. If included, distinguish block-level
   calibration from the strict equal-per-component definition of
   commensurability in Section 3.
3. **MDS instructions (Section 5).** State that `mds = TRUE` is required to
   compute these diagnostics; `dims = 2` is the default dimension only when
   MDS is requested. The current explanation can be read as saying MDS is
   computed by default.
4. **Draft placeholders.** The clustering paragraph still contains
   “with some defaults?” and “ARI; reference”. Replace these before submission;
   no explicit Hubert–Arabie ARI reference was found in the current bibliography.
5. **Conclusion after removing the geometry illustration.** The sentence saying
   the penguins application showed how scale and association shape geometry
   should be aligned with what is now actually demonstrated.
6. **Computational details.** The present rendering environment is R 4.6.0,
   manydist 0.5.2, tidymodels 1.5.0, dplyr 1.2.1, and palmerpenguins 0.1.1.
   The manually written R/tidymodels/dplyr versions differ. Generate these
   values from the rendering session or update them at final submission.
7. **Minor cleanup.** There are grammar/formatting details such as “that do not
   suffer” for a singular distance, “the the”, and “algorithn(s)”. These
   were preserved rather than silently changing the author's revisions.
   Two blue-colour markup problems were repaired after visual inspection:
   literal equation-label text and stray braces beside the association-aware
   equation reference.

## Completed in this task

- Removed `spec_type` from specification grids and benchmark output. Named
  presets select predefined settings; `preset = "custom"` selects explicit
  method settings. A legacy `spec_type` input column is ignored.
- Removed the redundant catalogue row for the generic `custom` entry. The
  full grid now contains 482 candidates, including the explicit custom
  combinations. Named-preset-only grids do not include that generic entry.
- Updated the source vignettes, help pages, tests, and the paper's benchmark
  paragraph, including the function spelling `all_dist_method_specs()` and a
  dynamically computed candidate count.
- All five package test files pass; the updated package was installed locally.
- The shortened 34-page manuscript rebuilt successfully. Pages 8 and 18–19
  were visually checked, and no literal unresolved equation labels remain.
- Comparing the source with the pre-edit snapshot confirmed preservation of
  the author's other current article edits.

The distance-algorithm issues above remain open and are not covered by the
passing API-regression tests. No Michel-comment status was changed.
