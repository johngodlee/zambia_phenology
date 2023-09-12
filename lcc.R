# Process Zambia land cover data
# John Godlee (johngodlee@gmail.com)
# 2021-07-23

# Packages
library(terra)

# Find files
lcc_files <- list.files("./dat/zambia_landcover", "*.hdf", full.names = TRUE)

# Construct layer definition
lcc_layers <- paste0("HDF4_EOS:EOS_GRID:", lcc_files, 
  ":MCD12Q1:LC_Type1")

# Mosaic
lcc_layers_list <- lapply(lcc_layers, rast)

lcc_mosaic <- mosaic(sprc(lcc_layers_list))

# Write to tif file
writeRaster(lcc_mosaic, "./dat/zambia_landcover/lcc.tif", overwrite = TRUE)
