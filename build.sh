#!/usr/bin/env sh

# Each script saves a .rds file which can then be loaded by the next script

# Extract plot data from SEOSAW database
Rscript data_prep.R  # 764

# Get raw values from various remote-sensing data
#Rscript modis_get.R
#Rscript lai_get.R
#Rscript trmm_get.R
#Rscript vipphen.R

# Extract statistics from raw remote-sensing data
#Rscript lai_extract.R
Rscript modis_extract.R  # 763
Rscript trmm_extract.R  # 762

# Calculate diversity metrics
Rscript diversity.R  # 710

# Run models
#Rscript analysis.R

