# JSS paper source and replication

This directory is the repository source of truth for the JSS manuscript.
The Dropbox/Overleaf project is treated as a publication mirror.

## Ownership

- Package API facts come from roxygen comments in ../../manydist/R.
- Method and preset identifiers come from manydist::dist_methods_tbl().
- Paper-only method formulas live in data/method-metadata.csv.
- The survey of related software lives in data/related-software.csv.
- article.qmd owns the manuscript narrative.
- replicate.R owns generated numerical results, tables, figures, and the
  replication manifest.
- dist/overleaf is generated and must not be edited by hand.

## First-time workflow

From the repository root:

~~~sh
Rscript paper/jss/validate.R
Rscript paper/jss/replicate.R
~~~

After installing the Quarto CLI:

~~~sh
quarto render paper/jss/article.qmd
Rscript paper/jss/bundle_overleaf.R
~~~

The optional _targets.R entry point runs the same replication function after
the R package "targets" has been installed. It is not required for the
zero-dependency replication entry point above.

## Generated files

The replication run writes only below generated/. Rendering writes below
dist/. Every generated artifact is included in generated/manifest.csv with
an MD5 checksum. The session information is stored alongside the numerical
results.

The first scaffold generates the shared penguin data, current distance
specifications, method registries, core distance summaries, and the MDS
figure. LOVO, KNN, and clustering targets can be migrated section by section
using MIGRATION.md.

