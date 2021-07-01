# Get MODIS 250 m (MOD13Q1) raw time series
# John Godlee (johngodlee@gmail.com)
# 2020-08-24

# Downloaded 4 granules
# 2010 - present
# 23 scenes per year per granule, except 2000 and 2021
# 1936 files total
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

# Import plot data 
plots <- readRDS("dat/sites_loc.rds")

plots_sf <- st_as_sf(plots, 
  coords = c("longitude_of_centre", "latitude_of_centre"), 
  crs = 4326)

# Get HDF file names 
files <- list.files("/Volumes/TOSHIBA_EXT/modis_250", pattern = "*.hdf", 
  full.names = TRUE)

# extract granule names
granExtract <- function(x) { 
  str_extract(x, "h[0-9][0-9]v[0-9][0-9]")
}

files_split <- split(files, granExtract(files))

# Create directory for output
data_dir <- "dat"
out_dir <- file.path("/Volumes/TOSHIBA_EXT/modis_250/evi_tif")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Define parameters 
nodata <- -3000
scaledata <- 0.0001

# For each granule create tif files
for (i in seq(length(files_split))) {

  gran_dir <- file.path(out_dir, names(files_split[i]))
  if (!dir.exists(gran_dir)) {
    dir.create(gran_dir)
  }

  # Get band names
  evi_list <- paste0("HDF4_EOS:EOS_GRID:", files_split[[i]], 
    ":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI")

  # Get .tif names
  out_names <- gsub("\\.hdf$", ".tif", 
    gsub(":.*", "", basename(evi_list)))

  # For each scene within a granule, create tif files 
  for (j in seq(length(evi_list))) {
    out_name <- file.path(gran_dir, out_names[j])
    print(sprintf("%s/%s : %s/%s : %s", 
        i, length(files_split), 
        j, length(evi_list), 
        evi_list[j]))
    if (!file.exists(out_name)) {
      gdal_translate(evi_list[j], 
        dst_dataset = out_name, 
        a_nodata = nodata)
    }
  }
}

# For each granule extract plot values
for (i in seq(length(files_split))) {
  # Define directory holding tif files for this granule
  gran_dir <- file.path(out_dir, names(files_split[i]))
  if (!dir.exists(gran_dir)) {
    dir.create(gran_dir)
  }

  # Read in .tif files
  tif_list <- lapply(list.files(gran_dir, "*.tif", full.names = TRUE), raster)

  # Add no-data field
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
  plots_fil <- st_filter(plots_sf, tif_poly)

  # For each scene within a granule, extract plot values 
  for (j in seq(length(tif_list))) {
    out <- as.character(raster::extract(tif_list[[j]], plots_fil, method = "bilinear"))
    outfile <- file.path(gran_dir, 
      paste0(names(tif_list)[j], "_", granExtract(files_split[[i]][1]), ".txt"))

    write_lines(out, outfile)
  }

  # Write file with plot cluster IDs
  plot_cls_file <- file.path(gran_dir, paste0(granExtract(files_split[[i]][1]), "_plot_id.txt"))
  write_lines(plots_fil$plot_cluster, plot_cls_file)

  # Remove .tif files
  #file.remove(file.path(out_dir, out_names))

  # Import .txt files, make columns in dataframe
  extract_files <- list.files(gran_dir, "*.txt", full.names = TRUE)
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
  write.csv(extract_df_clean, file.path(data_dir, 
      paste0("evi_", granExtract(files_split[[i]][1]), "_extract.csv")), 
    row.names = FALSE)
    
  # Remove .txt files
  #file.remove(list.files(out_dir, "*.txt", full.names = TRUE))
}

# Read .csv files per granule 
csv_files <- list.files(data_dir, "evi_h.*_extract.csv", full.names = TRUE)
csv_list <- lapply(csv_files, read.csv)

# Check that no plots have been included in both granules 
stopifnot(!any(
  intersect(unique(csv_list[[1]]$plot_cluster), unique(csv_list[[2]]$plot_cluster)),
  intersect(unique(csv_list[[1]]$plot_cluster), unique(csv_list[[3]]$plot_cluster)),
  intersect(unique(csv_list[[1]]$plot_cluster), unique(csv_list[[4]]$plot_cluster)),
  intersect(unique(csv_list[[2]]$plot_cluster), unique(csv_list[[3]]$plot_cluster)),
  intersect(unique(csv_list[[2]]$plot_cluster), unique(csv_list[[4]]$plot_cluster)),
  intersect(unique(csv_list[[3]]$plot_cluster), unique(csv_list[[4]]$plot_cluster))
  ))

# Join csv files
evi_ts_df <- do.call(rbind, csv_list)
evi_ts_df$date <- unlist(lapply(evi_ts_df$date, as.Date, 
    tryFormats = c("%Y-%m-%d", "%d/%m/%Y")))

# Write large raw data
saveRDS(evi_ts_df, "dat/evi.rds")

# Remove .csv files
#file.remove(list.files(data_dir, "evi_h.*_extract.csv", full.names = TRUE))


