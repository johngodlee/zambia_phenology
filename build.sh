#!/usr/bin/env sh

# Each script saves a .rds file which can then be loaded by the next script

Rscript data_prep.R
Rscript modis_get.R
Rscript modis_extract.R
Rscript diversity.R
Rscript vipphen.R
Rscript analysis.R

sed -i '7s/.*/Response \& DoF \& $\\delta$AIC\\textsubscript{m} \& $\\delta$AIC\\textsubscript{c} \& $\\chi$LRT \& p-LRT \\\\/g' out/spamm_stat.tex

sed -i 's/dT/$\\delta$T/g' out/hier_part.tex
