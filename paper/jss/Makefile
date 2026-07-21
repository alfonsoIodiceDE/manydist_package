.PHONY: validate replicate render bundle

validate:
	Rscript validate.R

replicate:
	Rscript replicate.R

render: replicate
	quarto render article.qmd

bundle: render
	Rscript bundle_overleaf.R

