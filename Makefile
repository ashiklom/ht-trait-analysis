.PHONY: all targets pdf

all: pdf

targets:
	Rscript -e "targets::tar_make()"

pdf:
	latexmk -cd -pdf manuscript/manuscript.tex
