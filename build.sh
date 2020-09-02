#!/usr/bin/env sh

Rscript data_prep.R
Rscript diversity.R
Rscript vipphen.R
Rscript modis_get.R
Rscript modis_extract.R
Rscript analysis.R
