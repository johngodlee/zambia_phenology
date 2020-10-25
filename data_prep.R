# Filter SEOSAW plot data for phenology analysis
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(dplyr)
library(tidyr)
library(sf)
library(tibble)

source("functions.R")

# Import data
plots <- read.csv("~/git_proj/seosaw_data/data_out/plots_v2.7.csv")

# Subset plots and data fields - Zambia ILUAii data only
plots_fil_sf <- plots %>% 
  filter(prinv == "Siampale A.") %>%  # Filter ILUAii plots
  dplyr::select(
    plot_id, plot_cluster, longitude_of_centre, latitude_of_centre, 
    richness, shannon, simpson, evenness, n_stems_gt10_ha, ba_ha, agb_ha,
    mat = bio1, diurnal_temp_range = bio2, map = bio12,
    clay = CLYPPT, sand = SNDPPT, cec = CECSOL,
    last_census_date) %>%  # Select columns
  group_by(plot_cluster) %>%  # Group by plot cluster
  summarise(
    plot_id = paste0(plot_id, collapse = ","),
    longitude_of_centre = mean(longitude_of_centre, na.rm = TRUE),
    latitude_of_centre = mean(latitude_of_centre, na.rm = TRUE),
    n_stems_gt10_ha = mean(n_stems_gt10_ha, na.rm = TRUE),
    ba_ha = mean(ba_ha, na.rm = TRUE),
    agb_ha = mean(agb_ha, na.rm = TRUE),
    mat = mean(mat, na.rm = TRUE),
    diurnal_temp_range = mean(diurnal_temp_range, na.rm = TRUE),
    map = mean(map, na.rm = TRUE),
    clay = mean(clay, na.rm = TRUE),
    sand = mean(sand, na.rm = TRUE),
    cec = mean(cec, na.rm = TRUE),
    last_census_date = first(last_census_date)) %>%  # Summarise by plot cluster
  st_as_sf(., coords = c("longitude_of_centre", "latitude_of_centre")) %>%  # Make spatial
  `st_crs<-`(4326) %>%  # Assign CRS
  mutate(plot_id_vec = strsplit(as.character(plot_id), split = ",")) %>%  # Make column of plot ID vectors 
  mutate(plot_id_length = sapply(.$plot_id_vec, length))  # Get number of plots in cluster

# Get statistics about sites
census <- unique(plots_fil_sf$last_census_date) 

n_total_sites <- plots %>%
  filter(prinv == "Siampale A.") %>%
  pull(plot_cluster) %>%
  unique() %>%
  length()

plots_fil_sf <- plots_fil_sf %>% 
  dplyr::select(-last_census_date)

# Create plot ID / plot Cluster lookup table
plot_id_lookup <- plots %>% 
  filter(plot_id %in% unlist(plots_fil_sf$plot_id_vec)) %>%
  dplyr::select(plot_cluster, plot_id)

# Write files
saveRDS(plots_fil_sf, "dat/plots.rds")
saveRDS(plot_id_lookup, "dat/plot_id_lookup.rds")

# Write stats to .tex
write(
  c(commandOutput(census, "censusDate"),
    commandOutput(n_total_sites, "nTotalSites")),
  file = "out/data_prep_vars.tex")
