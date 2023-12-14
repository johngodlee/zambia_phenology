#!/usr/bin/env sh

# Diff two drafts of the manuscript
# $1 = Original
# $2 = Revised version
# $3 = output (diff.tex)

latexdiff $1 $2 > $3
