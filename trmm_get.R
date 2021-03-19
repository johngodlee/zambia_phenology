# TRMM data processing
# John Godlee (johngodlee@gmail.com)
# 2020-10-13

# Set working directory

# Packages
library(sf)
library(dplyr)
library(raster)
library(ncdf4)
library(stringr)
library(readr)
library(tidyr)

# Import plot data
dat <- readRDS("dat/plots.rds")

# Get HDF file names 
files <- list.files("/Volumes/TOSHIBA EXT/trmm", pattern = "*.nc4$", 
  full.names = TRUE)

# Create directory for output
data_dir <- "dat"
out_dir <- file.path(data_dir, "trmm_tif")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Get .tif names
out_names <- gsub("\\.nc4$", ".tif", basename(files))

# Define bounding box to crop to
zam_bbox <- st_bbox(dat)

# For each scene, create tif files 
for (i in seq_along(files)) {
  out_name <- file.path(out_dir, out_names[i])

  print(sprintf("%s/%s : %s", i, length(files), out_name))

  if (!file.exists(out_name)) {
    nc_file <- nc_open(files[i])
    precip_array <- ncvar_get(nc_file, "precipitationCal")
    fillvalue <- ncatt_get(nc_file, "precipitationCal", "_FillValue")
    lon <- ncvar_get(nc_file, "lon")
    lat <- ncvar_get(nc_file, "lat")
    nc_close(nc_file) 
    precip_array[precip_array == fillvalue$value] <- NA

    r <- raster(precip_array, 
      xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
      crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    
    r <- flip(r, direction='y')

    r_crop <- crop(r, zam_bbox[c(1,3,2,4)])

    writeRaster(r_crop, out_name, overwite = TRUE)
  }
}

# Read in .tif files
tif_list <- lapply(list.files(out_dir, "*.tif", full.names = TRUE), raster)

# Rename stack layers to acquisition date
names(tif_list) <- unlist(lapply(tif_list, function(x) {
    str_extract(names(x), "[0-9]{8}")
}))

# For each scene, extract plot values 
lapply(seq(length(tif_list)), function(j) {
  print(sprintf("%s/%s : %s", j, length(tif_list), names(tif_list)[j]))
  out <- as.character(raster::extract(tif_list[[j]], dat, method = "bilinear"))
  outfile <- file.path(out_dir, "ext", paste0(names(tif_list)[j], ".txt"))

  write_lines(out, outfile)
})

# Write file with plot cluster IDs
plot_cls_file <- file.path(out_dir, "plot_id.txt")
write_lines(dat$plot_cluster, plot_cls_file)

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
saveRDS(extract_df_clean, file.path(data_dir, "trmm.rds"))
  
# Remove .txt files
file.remove(list.files(out_dir, "*.txt", full.names = TRUE))

# Remove .tif files
file.remove(file.path(out_dir, out_names))



