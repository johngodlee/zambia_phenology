# Get MODIS 250 m (MOD13Q1) raw time series
# John Godlee (johngodlee@gmail.com)
# 2020-08-24

# Downloaded 4 granules
# 2015 - present
# 22 scenes per year per granule.
# 517 files total
# MOD13Q1.A2015001.h20v10.006.2015295100845.hdf 
# PRODUCT.DATEAQUI.GRANUL.VER.DATEOFPRODUCT.hdf

# Set working directory

# Packages
library(rgdal)
library(stringr)
library(gdalUtils)
library(sf)
library(tidyr)
library(dplyr)
library(raster)
library(readr)

source("functions.R")

# Import plot data 
plots <- readRDS("dat/plots.rds")

# Get HDF file names 
files <- list.files("/Volumes/UNTITLED/modis_250", pattern = "*.hdf", 
  full.names = TRUE)

# Split file names by granule
granExtract <- function(x) { 
  str_extract(x, "h[0-9][0-9]v[0-9][0-9]")
}

files_split <- split(files, granExtract(files))

# Define CRS
crs_sinus <- CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

# Create directory for output
data_dir <- "dat"
out_dir <- file.path(data_dir, "tif")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Define parameters 
nodata <- -3000
scaledata <- 0.000001

# For each granule:
for (i in seq(length(files_split))) {

  # Get band names
  band_list <- lapply(files_split[[i]], get_subdatasets)
  evi_list <- unlist(lapply(band_list, `[`, 2))

  # Get .tif names
  out_names <- gsub("\\.hdf$", ".tif", 
    gsub(":.*", "", basename(evi_list)))

  # For each scene within a granule, create tif files 
  for (j in seq(length(evi_list))) {
    gdal_translate(evi_list[j], 
      dst_dataset = file.path(out_dir, out_names[j]), 
      a_nodata = nodata)
    }

  # Read in .tif files
  tif_list <- lapply(list.files(out_dir, "*.tif", full.names = TRUE), raster)
  lapply(tif_list, function(x) {
    NAvalue(x) <- nodata
      })

  # Rename stack layers to acquisition date
  names(tif_list) <- unlist(lapply(tif_list, function(x) {
      str_extract(names(x), "A[0-9]{7}")
  }))

  # Filter plots by granule coverage
  tif_poly <- st_as_sfc(st_bbox(tif_list[[1]])) %>%
    st_transform(., 4326)
  plots_fil <- st_filter(plots, tif_poly)

  # For each scene within a granule, extract plot values 
  for (j in seq(length(tif_list))) {
    out <- as.character(raster::extract(tif_list[[j]], plots_fil, method = "bilinear"))
    outfile <- file.path(out_dir, paste0(names(tif_list)[j], "_", granExtract(files_split[[i]][1]), ".txt"))

    write_lines(out, outfile)
  }

  # Write file with plot cluster IDs
  plot_cls_file <- file.path(out_dir, paste0(granExtract(files_split[[i]][1]), "_plot_id.txt"))
  write_lines(plots_fil$plot_cluster, plot_cls_file)

  # Remove .tif files
  file.remove(file.path(out_dir, out_names))

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
    gather(scene, evi, -plot_cluster, -granule) %>%
    mutate(date = as.Date(scene, "A%Y%j"),
      evi = as.numeric(evi) * scaledata) 

  # Write .csv for granule
  write.csv(extract_df_clean, file.path(data_dir, paste0(granExtract(files_split[[i]][1]), "_extract.csv")), row.names = FALSE)
    
  # Remove .txt files
  file.remove(list.files(out_dir, "*.txt", full.names = TRUE))
}

