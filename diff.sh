#!/usr/bin/env sh

# Diff two drafts of the manuscript
latexdiff phenology_old.tex phenology.tex > diff.tex

latexmk diff.tex
