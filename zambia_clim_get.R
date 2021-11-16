# Get Zambia climatic information
# John Godlee (johngodlee@gmail.com)
# 2021-07-22

# Packages
library(sf)
library(scico)
library(dplyr)
library(raster)
library(ggplot2)
library(gridExtra)

# Import data
bioclim <- stack("/Volumes/seosaw_spat/wc2.1_30s_bio/bioclim.vrt")

af <- st_read("dat/africa_countries/africa.shp")

# Extract Zambia from Africa
zambia <- af %>% 
  filter(sov_a3 == "ZMB")

# Get relevant layers from bioclim - 1, 12
bioclim_fil <- bioclim[[c(1,12)]]

# Clip to Zambia
bioclim_zambia <- mask(crop(bioclim_fil, zambia), zambia)

# Write to file
saveRDS(bioclim_zambia, "dat/bioclim_zambia.rds")

# Plot of MAT and MAP, for reference
pdf(file = "img/zambia_clim.pdf", width = 10, height = 15)
par(mfrow = c(2,1))
plot(bioclim_zambia$bioclim.1, main = "MAT")
plot(bioclim_zambia$bioclim.12, main = "MAP")
dev.off()


