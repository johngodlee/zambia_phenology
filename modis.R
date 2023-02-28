# Extract MODIS land cover dynamics 500 m (MCD12Q2) metrics for each site
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-02-24

# MODIS/Terra+Aqua Land Cover Dynamics Yearly L3 Global 500 m SIN Grid
# Downloaded 4 granules
# 2001-2021
# Files are yearly
# 88 files total
# MCD12Q2.A2021001.h21v10.061.2022294172721.hdf
# PRODUCT.DATEAQUI.GRANUL.VER.DATEOFPRODUCT.hdf

##' SDS Name - Description - Units - Data Type - Fill Value - No Data Value - Valid Range - Scale Factor
##' NumCycles - Total number of valid vegetation cycles with peak in product year - Number - 16-bit signed integer - 32767 - N/A - 1 to 7 - N/A
##' Greenup¹ - Date when EVI2 first crossed 15% of the segment EVI2 amplitude - Day - 16-bit signed integer - 32767 - N/A - 11138 to 32766 - N/A
##' MidGreenup¹ - Date when EVI2 first crossed 50% of the segment EVI2 amplitude - Day - 16-bit signed integer - 32767 - N/A - 11138 to 32766 - N/A
##' Peak¹ - Date when EVI2 reached the segment maximum - Day - 16-bit signed integer - 32767 - N/A - 11138 to 32766 - N/A
##' Maturity¹ - Date when EVI2 first crossed 90% of the segment EVI2 amplitude - Day - 16-bit signed integer - 32767 - N/A - 11138 to 32766 - N/A
##' Senescence¹ - Date when EVI2 last crossed 90% of the segment EVI2 amplitude - Day - 16-bit signed integer - 32767 - N/A - 11138 to 32766 - N/A
##' MidGreendown¹ - Date when EVI2 last crossed 50% of the segment EVI2 amplitude - Day - 16-bit signed integer - 32767 - N/A - 11138 to 32766 - N/A
##' Dormancy¹ - Date when EVI2 last crossed 15% of the segment EVI2 amplitude - Day - 16-bit signed integer - 32767 - N/A - 11138 to 32766 - N/A
##' EVI_Minimum - Segment minimum EVI2 value - EVI2 - 16-bit signed integer - 32767 - N/A - 0 to 10000 - 0.0001
##' EVI_Amplitude - Segment maximum - minimum EVI2 - EVI2 - 16-bit signed integer - 32767 - N/A - 0 to 10000 - 0.0001
##' EVI_Area - Sum of daily interpolated EVI2 from Greenup to Dormancy - EVI2 - 16-bit signed integer - 32767 - N/A - 0 to 3700 - 0.1
##' QA_Overall - QA code for entire segment - Class - 16-bit signed integer - 32767 - N/A - 0 to 3 - N/A
##' QA_Detailed - Bit-packed, SDS-specific QA codes - Bit Field - 16-bit signed integer - 32767 - N/A - 0 to 16383 - N/A

# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(terra)
library(sf)
library(parallel)

#' Fast rbind for dataframes
#'
#' @param x list of dataframe objects
#'
#' @return dataframe
#' 
fastRbind <- function(x) { 
  list2DF(lapply(setNames( seq_along(x[[1]]), names(x[[1]])), 
      function(i) {
        unlist(lapply(x, `[[`, i), FALSE, FALSE)
      }))
}

# Import plot data 
plots <- readRDS("dat/plots.rds")
trmm <- readRDS("dat/trmm.rds")
af <- st_read("dat/africa_countries/africa.shp")

# Get MODIS HDF file names 
files <- list.files("~/Desktop/MCD12Q2_2023-02-23", 
  pattern = "*.hdf", full.names = TRUE)

stopifnot(length(files) > 0)

granExtract <- function(x) { 
  str_extract(x, "h[0-9][0-9]v[0-9][0-9]")
}

# Split files list by granule
files_split <- split(files, granExtract(files))

# Define band names to extract 
band_names <- c("Greenup", "Maturity", "Senescence", "Dormancy", 
  "EVI_Area", "QA_Overall")

# For each granule extract band values for each plot
out <- fastRbind(lapply(seq_along(files_split), function(x) {
  # Read in first band to subset plots
  ex_band <- rast(files_split[[x]][1])

  # Filter plots by granule coverage
  tif_poly <- ex_band %>% 
    st_bbox() %>% 
    st_as_sfc()

  plots_fil <- plots %>% 
    st_transform(., crs(ex_band)) %>% 
    st_filter(., tif_poly)

  # For each band
  fastRbind(lapply(band_names, function(y) {
    # Generate band names
    band_list <- paste0("HDF4_EOS:EOS_GRID:", files_split[[x]], 
      ":MCD12Q2:", y)

    # Read in rasters 
    rast_list <- lapply(band_list, function(z) {
      rast(z)[[1]]
    })

    # Rename stack layers to acquisition date
    names(rast_list) <- gsub("A", "", str_extract(band_list, "A[0-9]{4}"))

    fastRbind(mclapply(seq_along(rast_list), function(z) {
      message(names(files_split)[x], ":", y, ":", names(rast_list)[z])

      ext_val <- unname(terra::extract(rast_list[[z]], st_coordinates(plots_fil)))

      data.frame(
        plot_cluster = plots_fil$plot_cluster,
        var = y,
        year = names(rast_list[z]),
        ext_val)
    }, mc.cores = detectCores()-1))
  }))
}))

# Remove duplicate entries by averaging their values
out_clean <- out %>% 
  group_by(plot_cluster, year, var) %>% 
  summarise(ext_val = mean(ext_val))

# Check that all plots have correct data extracted
stopifnot(all(
    table(paste(out_clean$plot_cluster, out_clean$year, sep = "_")) == 
      length(band_names)))

# Spread dataframe
out_spread <- out_clean %>% 
  pivot_wider(
    names_from = var,
    values_from = ext_val)

# Calculate length of greening and senescence period
out_spread$green_rate <- out_spread$Maturity - out_spread$Greenup
out_spread$senes_rate <- out_spread$Dormancy - out_spread$Senescence
out_spread$season_length <- out_spread$Dormancy - out_spread$Greenup

# Check all scenese of "best" quality
stopifnot(all(out_spread$QA_Overall == 0 | is.na(out_spread$QA_Overall)))

# Calculate average values and lag times
out_avg <- out_spread %>% 
  mutate(
    across(all_of(band_names[1:4]), 
      ~as.Date(.x, origin = "1970-01-01") - as.Date(paste0(year, "-01-01"))),
    across(all_of(band_names[1:4]), 
      ~ifelse(.x > 182, .x-365, .x))) %>% 
  left_join(., trmm, by = c("plot_cluster", "year")) %>%
  mutate(
    start_lag = -(Greenup - trmm_start),
    end_lag = Dormancy - trmm_end) %>% 
  ungroup() %>% 
  dplyr::select(-QA_Overall, -year) %>%
  group_by(plot_cluster) %>% 
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

# Check dataframe contains all plots and no duplicates
stopifnot(all(out_avg$plot_cluster %in% plots$plot_cluster))
stopifnot(all(!duplicated(out_avg$plot_cluster)))

# Write extracted values to file
saveRDS(out_avg, "dat/modis.rds")
