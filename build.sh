#!/usr/bin/env sh

# Each script saves a .rds file which can then be loaded by the next script

# Extract plot data from SEOSAW database
Rscript data_prep.R

# Get raw values from various remote-sensing data
#Rscript modis_get.R
#Rscript lai_get.R
#Rscript trmm_get.R

# Extract statistics from raw remote-sensing data
Rscript modis_extract.R
Rscript vipphen.R
#Rscript lai_extract.R
Rscript trmm_extract.R

# Calculate diversity metrics
Rscript diversity.R

# Run models
Rscript analysis.R

#sed -i '7s/.*/Response \& DoF \& $\\delta$AIC\\textsubscript{m} \& $\\delta$AIC\\textsubscript{c} \& $\\chi$LRT \& p-LRT \\\\/g' out/spamm_stat.tex

#sed -i 's/dT/$\\delta$T/g' out/hier_part.tex
