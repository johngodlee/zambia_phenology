# Extract MODIS VIPPHEN EVI metrics
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)

# Import data
plots <- readRDS("dat/plots.rds")

af <- st_read("dat/africa_countries/africa.shp")

zambia <- af %>% 
  filter(sov_a3 == "ZMB")

# Create plot extent
plots_extent <- extent(plots) 
plots_extent[c(1,3)] <- plots_extent[c(1,3)] -1
plots_extent[c(2,4)] <- plots_extent[c(2,4)] +1

plots_bbox <- st_as_sfc(st_bbox(plots_extent)) %>%
  `st_crs<-`(4326)

# List of vipphen layers
vipphen_file_list <- list.files(path = "/Volumes/UNTITLED/modis_vipphen", 
  pattern = ".*tif$", full.names = TRUE)

# Split by variable
var_list <- as.numeric(gsub("\\.tif", "", gsub(".*_", "", vipphen_file_list)))
vipphen_file_list_split <- split(vipphen_file_list, var_list)

# Stack by variable
phen_stack_list <- lapply(vipphen_file_list_split, stack)

# Select variables I care about
phen_stack_list_fil <- phen_stack_list[c(1,2,3,5,6,22,23,24,25)]

# Crop phenology
phen_crop_list <- lapply(phen_stack_list_fil, crop, y = plots_extent)

# Set CRS of raster to WGS84
phen_wgs_list <- lapply(phen_crop_list, projectRaster, 
  crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

phen_zam_list <- lapply(phen_wgs_list, mask, zambia)
phen_zam_mean_stack <- stack(lapply(phen_zam_list, calc, mean))

saveRDS(phen_zam_mean_stack, "dat/vipphen_stack.rds")
saveRDS(phen_stack_list[[26]], "dat/phen_zam_reliab.rds")

# Pull values and take mean across years
plots_phen_extract <- as.data.frame(do.call(cbind, lapply(phen_wgs_list, function(x) { 
      raster::extract(x, plots, method = "simple") %>%
  as.data.frame(.) %>%
  mutate_all(~ case_when(
      . == -1 ~ NA_real_,
      . < -5000 ~ NA_real_,
      TRUE ~ .)) %>%
  rowMeans(., na.rm = TRUE)
})))

names(plots_phen_extract) <- paste0("vipphen_", c("s1_start", "s1_end", "s1_length", 
  "s1_green_rate", "s1_senes_rate", "cum_vi", "avg_vi", "bg_vi", "n_seasons"))

plots_vipphen <- plots %>%
  dplyr::select(plot_cluster) %>%
  st_drop_geometry() %>%
  bind_cols(., plots_phen_extract) 

saveRDS(plots_vipphen, "dat/vipphen.rds")

