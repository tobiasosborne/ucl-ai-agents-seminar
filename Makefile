# Makefile for UCL AI Agents Seminar slides
# Requires: lualatex (LuaHBTeX 1.17+)

LATEX    := lualatex
LATEXOPT := --interaction=nonstopmode --halt-on-error
SRCDIR   := slides

SECTIONS := section00 section01 section02 section03 section04
PDFS     := $(SECTIONS:%=$(SRCDIR)/%.pdf)

# Theme files (rebuild if any change)
THEME := $(wildcard $(SRCDIR)/beamer*TJO.sty) $(SRCDIR)/preamble.tex

# Destination directories
DEST_00 := 00-why/slides.pdf
DEST_01 := 01-stateless-nondeterministic-function/slides.pdf
DEST_02 := 02-illusion-of-chat/slides.pdf
DEST_03 := 03-primitive-agent/slides.pdf
DEST_04 := 04-audience-vote/slides.pdf

.PHONY: all clean dist seminar s00 s01 s02 s03 s04

all: seminar

# Build the unified deck
seminar: $(SRCDIR)/seminar.pdf
	@echo "Unified deck built: $(SRCDIR)/seminar.pdf"

# Compile each section (two passes for correct slide numbering)
$(SRCDIR)/%.pdf: $(SRCDIR)/%.tex $(THEME)
	cd $(SRCDIR) && $(LATEX) $(LATEXOPT) $(<F)
	cd $(SRCDIR) && $(LATEX) $(LATEXOPT) $(<F)

# Copy PDFs to section directories
dist: $(PDFS)
	cp $(SRCDIR)/section00.pdf $(DEST_00)
	cp $(SRCDIR)/section01.pdf $(DEST_01)
	cp $(SRCDIR)/section02.pdf $(DEST_02)
	cp $(SRCDIR)/section03.pdf $(DEST_03)
	cp $(SRCDIR)/section04.pdf $(DEST_04)
	@echo "All slides built and distributed."

clean:
	cd $(SRCDIR) && rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.pdf

# Build individual sections
s00: $(SRCDIR)/section00.pdf
	cp $(SRCDIR)/section00.pdf $(DEST_00)
s01: $(SRCDIR)/section01.pdf
	cp $(SRCDIR)/section01.pdf $(DEST_01)
s02: $(SRCDIR)/section02.pdf
	cp $(SRCDIR)/section02.pdf $(DEST_02)
s03: $(SRCDIR)/section03.pdf
	cp $(SRCDIR)/section03.pdf $(DEST_03)
s04: $(SRCDIR)/section04.pdf
	cp $(SRCDIR)/section04.pdf $(DEST_04)
