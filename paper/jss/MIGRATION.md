# Migration checklist

The legacy manuscript reviewed for this scaffold is:

- Path: /Users/Alfo/Library/CloudStorage/Dropbox/Applicazioni/Overleaf/jss mixed distances/article.tex
- SHA-256 on 2026-07-21:
  34ebc14712be8c511964d4e39acb483170c7688a94d0a2fb067d948e59cdd0f1

Do not edit the legacy file during migration. Move one complete content block
at a time, regenerate its artifacts, compare the rendered result, and only
then mark it complete here.

## Status

- [ ] Front matter, abstract, and author addresses
- [x] Introduction and related-software tables
- [x] Unified mixed-variable distance framework
- [x] Categorical method formula table
- [ ] Package design and mdist examples
- [ ] Presets and method catalogue
- [ ] Benchmarking examples
- [ ] LOVO diagnostics and four LOVO figures
- [ ] step_mdist and KNN workflow
- [ ] PAM and spectral clustering workflow
- [ ] Spectral embedding figure
- [ ] MDS illustration
- [ ] Conclusion
- [ ] Bibliography and JSS class/template

## Migration rules

1. Use the current API names method_cat and method_num.
2. Use prepare_penguins() for the shared six-predictor data.
3. Use paper_distance_specifications() instead of repeating preset lists.
4. Never paste printed numerical output into article.qmd.
5. Every included table and figure must be created by run_replication().
6. Add an explicit seed for every stochastic computation.
7. Keep article.tex and the Overleaf bundle as generated products.
8. Prefer tidyverse functions in manuscript and replication code, including
   readr for delimited files and dplyr/purrr for data transformations.
