.PHONY: all targets

all: manuscript/manuscript.pdf

targets:
	Rscript -e "targets::tar_make()"

manuscript/manuscript.pdf:
	latexmk -cd -pdf manuscript/manuscript.tex
