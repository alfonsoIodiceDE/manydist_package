# Coauthor review progress

Updated 2026-10-27.

Progress copy: [Annotated PDF with colour-coded comments](article_cc_mvdv_progress.pdf).

The [original annotated PDF](article_cc_mvdv.pdf) is unchanged. A [backup of the progress copy before this update](review-history/article_cc_mvdv_progress_before-2026-10-27-8d5be1a2.pdf) preserves the author's manually added answer. The updated progress copy retains that answer too.

This review compares annotations with the saved `article.qmd` source. The annotated PDF still contains the original manuscript text and its original 36 pages; it is not a newly typeset version of the revised article. The manuscript PDF was not rebuilt in this colour-update step.

## Colour convention and current totals

- **Green: 26** — changes made in response to the comment, including removal from the visible paper even if helper code remains hidden.
- **Light blue: 1** — answered without changing the manuscript.
- **Original colours: 72** — five partial cases and 67 other pending comments.
- All 99 original annotation records remain editable, with reviewer names, comment bodies, rich text, and reply relationships preserved. The author's appended answer is retained verbatim. No annotations were deleted or flattened.
- Colour is an editorial progress assessment, not blanket approval of the surrounding passage. All previously queued colour updates have now been applied.

## Addressed — green

| Annotation | Reviewer | Original PDF page | Evidence |
|---|---|---:|---|
| R2-001 | Carlo | 1 | The author confirmed on 2026-10-27 that the abstract comment is addressed; its green colour has now been applied. |
| R2-003 | Carlo | 3 | The introduction now describes the progression from a direct distance-based learning interface to tidymodels workflows for clustering and classification. |
| R2-004 | Carlo | 5 | Both related-software CSVs now use cluspcamix, matching the existing prose. |
| R2-005 | Carlo | 8 | Table 4 includes u_dep_bw; Section 4.3 explains the numerical-block construction and gives its equivalent block_comm = TRUE settings. |
| R2-007 | Carlo | 11 | The presets introduction has been replaced by an explanation of specifying components versus selecting presets; the flagged wording has been removed. |
| R2-008 | Michel | 11 | The first distance example now explicitly states that species, stored in the first column, is excluded from the predictors. It uses the new default rather than retaining the struck unnormalised-Gower introduction. |
| R2-011 | Michel | 11 | The default and additional Gower examples now use penguin_train[, -1], as suggested, with the name-based alternative retained in a footnote. |
| R2-012 | Michel | 11 | The default remains the first worked example. At the author's request, a second, assignment-only mdist call now constructs d_gower_sum with range/matching/FALSE settings. The accompanying text distinguishes summed contributions from the variable-averaged Gower implementations and the gower preset, qualified to complete equally weighted data. No Gower print or equality-check output is displayed. |
| R2-013 | Michel | 11 | The flagged normalised-Gower comparison chunk, all.equal output, and following sentence have been removed. The separate MDist return-object explanation is now provided using d_default for R2-014. |
| R2-014 | Michel | 11 | The author chose to replace the Gower example with the current default. Added print(d_default) before Other distance specifications, explained the settings, and distinguished the MDist object from d_default$distance and its full matrix conversion. The exact example and full manuscript rendered successfully on 2026-10-27. |
| R2-015 | Michel | 11 | Replaced the Gower construction example with explicit d_default using method_num = none, method_cat = tvd, and logical commensurable = TRUE, followed by the equivalent argument-free call. Used the package's newly agreed actual default none rather than the earlier std suggestion. Executed and verified identical distances for both calls and u_mix, and rendered the updated manuscript. |
| R2-016 | Michel | 11 | The presets section now identifies preset = "gower" as the averaged Gower distance used by the referenced implementations. |
| R2-017 | Michel | 11 | The presets introduction now explicitly distinguishes component specifications from named presets, including non-additive constructions. |
| R2-019 | Carlo | 12 | The Table 4 caption now starts with the plural "Distances implemented via presets". |
| R2-020 | Michel | 12 | Removed the redundant d_u_mix \|> print() call from the presets example; the MDist object has already been illustrated with print(d_default). The author also confirmed the predictor selection should be [, -1], which has been corrected. |
| R2-021 | Carlo | 13 | The supervised-distance introduction now reads "In supervised settings". |
| R2-024 | Michel | 13 | The explanation names the fitted category dissimilarities and numerical scaling used for the training observations. It now matches the robust/tvd/FALSE example and explicitly states that no commensurability weights are estimated there. |
| R2-029 | Carlo | 14 | The additional-diagnostics example now uses "e.g.,". |
| R2-030 | Carlo | 14 | The same sentence now uses the plural "clusters". |
| R2-031 | Michel | 14 | The paragraph now begins "Additional diagnostics". |
| R2-032 | Michel | 14 | The mad_importance bullet now uses a colon rather than "is". |
| R2-033 | Michel | 14 | The relative_distance bullet now uses ": normalized version" rather than "is the normalized version". |
| R2-035 | Michel | 14 | The interpretive sentence no longer says "as the relative"; the requested "the" has been removed. |
| R2-036 | Michel | 14 | The interpretive sentence now uses "relative LOVO contributions" in the plural. |
| R2-045 | Michel | 16 | The lengthy formatting code is omitted from the visible manuscript through include: FALSE and echo: FALSE. This satisfies removal of the displayed code; physical deletion of a hidden helper chunk is not required. The broader package-default suggestion remains separate under R2-044. |
| R2-050 | Michel | 17 | The comparative LOVO distance plot now calls ggplot2::autoplot(lovo_comparison), using its default metric. |

## Replied without manuscript changes — light blue

### R2-006 — Michel — original PDF p. 9

The author replied in the PDF that the existing penguins data should not be duplicated merely to drop a variable and missing cases, and that this is a useful minimal-preparation example. Manuscript unchanged for this suggestion; mark light blue.

Author's answer appended to the PDF, verbatim:

> I would keep this, I would not define a dataset that is already availble just to drop one variable and some missings. Also, I believe it serves as an example of minimal data preparation

## Partial — original colour retained

- **R2-022** (p. 13): The response-aware training example now explicitly supplies method_cat = tvd and custom numerical/commensurability settings, with an association-aware qualification in the prose. The separate question about ignored responses and inconsistent preset warnings has not been resolved; retain the original colour.
- **R2-025** (p. 13): The training and new-data examples now consistently use robust/tvd/FALSE, with test responses excluded and matching prose. The broader request for a separate unsupervised new-data example and explicit matrix dimensions remains pending.
- **R2-044** (p. 16): The formatting code is hidden and the visible calls are shorter, but the suggested package-default/label changes have not been implemented or agreed.
- **R2-061** (p. 19): The package grid now uses none rather than std, but the manuscript's specification-grid paragraph still says std.
- **R2-074** (p. 22): The abstract was softened, but no leakage demonstration was restored and the overall author decision has not yet been explicitly confirmed.

All other comments remain open for the author's step-by-step review. [coauthor-review-status.json](coauthor-review-status.json) records the original annotation IDs, statuses, evidence, colour values, and source hashes.

## New checks and corrections in this pass

- Corrected the `u_mix` example to `penguin_train[, -1]`, as the author explicitly confirmed. No other manuscript changes were made in this pass.
- **Resolved after this colour update:** the author aligned both response-aware calls to `method_num = "robust", method_cat = "tvd", commensurable = FALSE`. The prose now refers to fitted category dissimilarities and robust numerical preprocessing, explicitly noting that no commensurability weights are estimated. The progress PDF colours are unchanged.
- The specification-grid paragraph still says `std` rather than `none`; R2-061 remains partial.
- The theoretical statement about preserving PCs' relative importance remains to be reconciled with whitening in the explicit blockwise construction. No mathematical wording was changed.

