# Makefile for GC Conference 2026 talk
# Requires: lualatex (LuaHBTeX 1.17+)

LATEX    := lualatex
LATEXOPT := --interaction=nonstopmode --halt-on-error

.PHONY: all clean

all: talk.pdf
	@echo "Talk built: talk.pdf"

talk.pdf: talk.tex preamble.tex $(wildcard beamer*TJO.sty)
	$(LATEX) $(LATEXOPT) talk.tex
	$(LATEX) $(LATEXOPT) talk.tex

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.pdf *.listing
