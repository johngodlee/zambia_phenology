# TRMM data processing
# John Godlee (johngodlee@gmail.com)
# 2020-10-13

# Set working directory

# Packages
library(sf)
library(dplyr)
library(raster)
library(gdalUtils)
library(stringr)
library(readr)
library(tidyr)

source("functions.R")

# Import plot data
dat <- readRDS("dat/plots.rds")

# Get HDF file names 
files <- list.files("/Volumes/UNTITLED/trmm", pattern = "*.nc4$", 
  full.names = TRUE)

# Create directory for output
data_dir <- "dat"
out_dir <- file.path(data_dir, "trmm_tif")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Get band names
band_list <- lapply(files, get_subdatasets)
trmm_list <- unlist(lapply(band_list, `[`, 1))

# Get .tif names
out_names <- gsub("\\.nc4$", ".tif", 
  gsub(":.*", "", basename(trmm_list)))

# Define bounding box to crop to
zam_bbox <- st_bbox(dat)

# For each scene, create tif files 
for (j in seq(length(trmm_list))) {
  gdal_translate(trmm_list[j], 
    dst_dataset = file.path(out_dir, out_names[j]))

  x <- raster(file.path(out_dir, out_names[j]))
  x <- t(flip(x, direction='y'))
  crs(x) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

  x <- crop(x, zam_bbox[c(1,3,2,4)])

  file.remove(file.path(out_dir, out_names[j]))

  writeRaster(x, file.path(out_dir, out_names[j]), overwite = TRUE)
}

# Read in .tif files
tif_list <- lapply(list.files(out_dir, "*.tif", full.names = TRUE), raster)

# Rename stack layers to acquisition date
names(tif_list) <- unlist(lapply(tif_list, function(x) {
    str_extract(names(x), "[0-9]{8}")
}))

# For each scene, extract plot values 
lapply(seq(length(tif_list)), function(j) {
  out <- as.character(raster::extract(tif_list[[j]], dat, method = "bilinear"))
  outfile <- file.path(out_dir, paste0(names(tif_list)[j], ".txt"))

  write_lines(out, outfile)
})

# Write file with plot cluster IDs
plot_cls_file <- file.path(out_dir, "plot_id.txt")
write_lines(dat$plot_cluster, plot_cls_file)

# Remove .tif files
file.remove(file.path(out_dir, out_names))

# Import .txt files, make columns in dataframe
extract_files <- list.files(out_dir, "*.txt", full.names = TRUE)
extract_vec <- lapply(extract_files, readLines)
extract_df <- as.data.frame(do.call(cbind, extract_vec))
names(extract_df) <- c(names(tif_list), "plot_cluster")
extract_df[,seq(length(extract_df) - 1)] <- lapply(
  extract_df[,seq(length(extract_df) - 1)], as.numeric)

# Make dataframe clean
extract_df_clean <- extract_df %>%
  gather(scene, precip, -plot_cluster) %>%
  mutate(date = as.Date(scene, "%Y%m%d"),
    precip = as.numeric(precip))

# Write .csv 
write.csv(extract_df_clean, file.path(data_dir, "trmm_extract.csv"), row.names = FALSE)
  
# Remove .txt files
file.remove(list.files(out_dir, "*.txt", full.names = TRUE))



