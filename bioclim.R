# Get Zambia long-term climate data
# John Godlee (johngodlee@gmail.com)
# 2021-07-22

# Packages
library(sf)
library(terra)

# Import bioclim data
bioclim <- rast(c(
    "./dat/wc2.1_30s_bio/wc2.1_30s_bio_1.tif",
    "./dat/wc2.1_30s_bio/wc2.1_30s_bio_12.tif"))

# Import plot data
plots <- readRDS("./dat/plots.rds")

# Import africa
af <- st_read("dat/africa_countries/africa.shp")

# Extract Zambia from Africa
zambia <- af[af$sov_a3 == "ZMB",]

# Convert Zambia to spatVector
zambia_sv <- vect(zambia)

# Clip to Zambia
bioclim_zambia <- mask(crop(bioclim, zambia_sv), zambia_sv)

# Pack for writing to file
bioclim_pack <- wrap(bioclim_zambia)

# Write to file
saveRDS(bioclim_pack, "dat/bioclim_zambia.rds")

# Plot of MAT and MAP, for reference
pdf(file = "img/zambia_clim.pdf", width = 10, height = 15)
par(mfrow = c(2,1))
plot(bioclim_zambia$wc2.1_30s_bio_1, main = "MAT")
plot(bioclim_zambia$wc2.1_30s_bio_12, main = "MAP")
dev.off()

# Add 10 km buffer around Zambia
zambia_buff <- st_buffer(zambia, 10000)

# Convert Zambia with buffer to spatVector
zambia_buff_sv <- vect(zambia_buff)

# Clip to Zambia
bioclim_zambia_buff <- mask(crop(bioclim, zambia_buff_sv), zambia_buff_sv)

# Extract BioClim values for each plot
bioclim_ext <- cbind(
  plots$plot_cluster, 
  terra::extract(bioclim_zambia_buff, st_coordinates(plots)))

names(bioclim_ext) <- c("plot_cluster", "mat", "map")

# Save to file
saveRDS(bioclim_ext, "./dat/bioclim.rds")

