# Process Zambia season length raster
# John L. Godlee (johngodlee@gmail.com)
# Last updated: 2023-09-17

# Packages
library(sf)
library(terra)

# Import data
af <- st_read("dat/africa_countries/africa.shp")
r <- rast("~/Downloads/modis_zambia.tif")
lcc <- rast("./dat/zambia_landcover/lcc.tif")

# Define pad value to get good dormancy value from mid-greendown
dorm_pad <- 30

# Find season length
slen <- r[[2]] - r[[1]]
values(slen) <- values(slen) + dorm_pad

# Project Zambia outline to same CRS as VRTs
zambia <- af[af$sov_a3 == "ZMB",]
zambia_trans <- st_transform(zambia, st_crs(slen))

# Mask season length layer by Zambia 
slen_mask <- mask(slen, vect(zambia_trans))

# Mask land cover to woodlands and savannas 
lcc_fil <- lcc
lcc_fil[lcc_fil %in% 1:10] <- NA
lcc_wgs84 <- project(lcc_fil, crs(zambia_trans))
zambia_sv <- vect(zambia_trans)
lcc_mask <- mask(crop(lcc_wgs84, zambia_sv), zambia_sv)

# Aggregate to same resolution as phenology
lcc_agg <- resample(lcc_mask, slen_mask)

# Mask phenology by land cover
slen_lcc <- mask(slen_mask, lcc_agg, inverse = TRUE)

# Write to file
writeRaster(slen_lcc, "./dat/slen_zambia.tif", overwrite = TRUE)

