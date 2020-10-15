# LAI extraction 
# John Godlee (johngodlee@gmail.com)
# 2020-10-09

# Set working directory

# Packages
library(dplyr)
library(tidyr)
library(sf)
library(raster)
library(stringr)
library(gdalUtils)
library(readr)

source("functions.R")

# Import plot data 
dat <- readRDS("dat/plots.rds")

# Get HDF file names 
files <- list.files("/Volumes/john/modis_lai", pattern = "*.hdf", 
  full.names = TRUE)

# Split by granule
files_split <- split(files, granExtract(files))

# Define sinusoidal CRS
crs_sinus <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# Define directory paths
data_dir <- "dat"
out_dir <- file.path(data_dir, "lai_tif")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

for (i in c(2,3,4)) {

  # Get band names
  lai_list <- paste0("HDF4_EOS:EOS_GRID:", files_split[[i]], 
    ":MOD_Grid_MOD15A2H:Lai_500m")

  # Get .tif names
  out_names <- gsub("\\.hdf$", ".tif", 
    gsub(":.*", "", basename(lai_list)))

  # For each scene within a granule, create tif files 
  for (j in seq(length(lai_list))) {
    gdal_translate(lai_list[j], 
      dst_dataset = file.path(out_dir, out_names[j]))
    }

  tif_list <- lapply(list.files(out_dir, "*.tif", full.names = TRUE), raster)

  names(tif_list) <- unlist(lapply(tif_list, function(x) {
      str_extract(names(x), "A[0-9]{7}")
  }))

  tif_poly <- st_as_sfc(st_bbox(tif_list[[1]])) %>%
    st_transform(., 4326)
  plots_fil <- st_filter(dat, tif_poly)

  for (j in seq(length(tif_list))) {
    out <- as.character(raster::extract(tif_list[[j]], plots_fil, method = "bilinear"))
    outfile <- file.path(out_dir, paste0(names(tif_list)[j], "_", granExtract(files_split[[i]][1]), ".txt"))

    write_lines(out, outfile)
  }

  # Write file with plot cluster IDs
  plot_cls_file <- file.path(out_dir, paste0(granExtract(files_split[[i]][1]), "_plot_id.txt"))
  write_lines(plots_fil$plot_cluster, plot_cls_file)

  # Import .txt files, make columns in dataframe
  extract_files <- list.files(out_dir, "*.txt", full.names = TRUE)
  extract_vec <- lapply(extract_files, readLines)
  extract_df <- as.data.frame(do.call(cbind, extract_vec))
  names(extract_df) <- c(names(tif_list), "plot_cluster")
  extract_df[,seq(length(extract_df) - 1)] <- lapply(
    extract_df[,seq(length(extract_df) - 1)], as.numeric)

  # Add granule ID
  extract_df$granule <- names(files_split)[i]

  # Make dataframe clean
  extract_df_clean <- extract_df %>%
    gather(scene, lai, -plot_cluster, -granule) %>%
    mutate(date = as.Date(scene, "A%Y%j"))

  # Write .csv for granule
  write.csv(extract_df_clean, file.path(data_dir, 
      paste0("lai_", granExtract(files_split[[i]][1]), "_extract.csv")), 
    row.names = FALSE)
    
  # Remove .txt files
  file.remove(list.files(out_dir, "*.txt", full.names = TRUE))

  # Remove .tif files
  file.remove(file.path(out_dir, out_names))
}

