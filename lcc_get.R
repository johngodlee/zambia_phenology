# Process Zambia land cover data
# John Godlee (johngodlee@gmail.com)
# 2021-07-23

# Packages
library(raster)
library(gdalUtils)

# Find files
lcc_files <- list.files("dat/zambia_landcover", "*.hdf", full.names = TRUE)

# Construct layer definition
lcc_layer_defs <- paste0("HDF4_EOS:EOS_GRID:", lcc_files, 
  ":MCD12Q1:LC_Type1")

# For each file, save layer as .tif
for (j in seq(length(lcc_files))) {
  out_name <- gsub("\\.hdf", ".tif", lcc_files[j])
  if (!file.exists(out_name)) {
    gdal_translate(lcc_layer_defs[j], 
      dst_dataset = out_name, 
      a_nodata = 255)
  }
}

# Import tifs
lcc_layers_files <- list.files("dat/zambia_landcover", "*.tif", full.names = TRUE)

# Mosaic
lcc_layers_list <- lapply(lcc_layers_files, raster)
lcc_layers_list$fun <- mean

lcc_mosaic <- do.call(mosaic, lcc_layers_list)

# Write to tif file
writeRaster(lcc_mosaic, "dat/zambia_landcover/lcc.tif")
