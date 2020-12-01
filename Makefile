.PHONY: all

all: manuscript/manuscript.pdf

manuscript/manuscript.pdf:
	latexmk -cd -pdf manuscript/manuscript.tex
