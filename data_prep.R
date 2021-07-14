# Filter SEOSAW plot data for phenology analysis
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Filter to ILUAii plots
# Summarise to site (x4 plots)
# Save site locations for modis_get.R and trmm_get.R
# Filter to stems >10 cm DBH, alive
# Summarise stems to trees
# Filter to sites with >5 species with >1 individual
# Filter to sites with >50 trees ha
# Filter to sites with no non-native species
# Create basal area abundance matrices

# Packages
library(dplyr)
library(tidyr)
library(sf)
library(tibble)

source("functions.R")

# Import data
plots <- read.csv("dat/plots_v2.12.csv")
plots_sp <- read.csv("dat/plots_spatial_v2.12.csv")

stems <- read.csv("dat/stems_iluaii_v2.12.csv")

# Subset plots and data fields - Zambia ILUAii data only
plots_iluaii <- plots %>% 
  filter(prinv == "Siampale A.")  # Filter to ILUAii plots

# Get statistics about sites
census <- unique(plots_iluaii$census_date) 

n_total_sites <- plots_iluaii %>%
  pull(plot_cluster) %>%
  unique() %>%
  length()

# Save site locations to file, for modis_get.R and trmm_get.R
plots_iluaii %>% 
  group_by(plot_cluster) %>%  # Group by plot cluster
  summarise(
    plot_id = paste0(plot_id, collapse = ","),
    longitude_of_centre = mean(longitude_of_centre, na.rm = TRUE),
    latitude_of_centre = mean(latitude_of_centre, na.rm = TRUE)) %>% 
  saveRDS(., "dat/sites_loc.rds")

# Summarise to sites
plots_fil_sf <- plots_iluaii %>%  
  left_join(., plots_sp, by = "plot_id") %>%
  dplyr::select(
    plot_id, plot_cluster, longitude_of_centre, latitude_of_centre, 
    plot_area, n_stems_ge10,
    richness, shannon, simpson, evenness, ba_ha, agb_ha,
    mat = bio1, diurnal_temp_range = bio2, map = bio12,
    clay = soil_clay, sand = soil_sand, cec = soil_cation_ex_cap, 
    census_date) %>%  # Select columns
  group_by(plot_cluster) %>%  # Group by plot cluster
  summarise(
    plot_id = paste0(plot_id, collapse = ","),
    longitude_of_centre = mean(longitude_of_centre, na.rm = TRUE),
    latitude_of_centre = mean(latitude_of_centre, na.rm = TRUE),
    census_date = first(census_date),
    ba_ha = mean(ba_ha, na.rm = TRUE),
    agb_ha = mean(agb_ha, na.rm = TRUE),
    n_stems_ge10_ha = mean(n_stems_ge10, na.rm = TRUE) / 
      sum(plot_area, na.rm = TRUE),
    mat = mean(mat, na.rm = TRUE),
    diurnal_temp_range = mean(diurnal_temp_range, na.rm = TRUE),
    map = mean(map, na.rm = TRUE),
    clay = mean(clay, na.rm = TRUE),
    sand = mean(sand, na.rm = TRUE),
    cec = mean(cec, na.rm = TRUE)) %>%  # Summarise by plot cluster
  st_as_sf(., coords = c("longitude_of_centre", "latitude_of_centre")) %>%  # Make spatial
  `st_crs<-`(4326) %>%  # Assign CRS
  mutate(plot_id_vec = strsplit(as.character(plot_id), split = ",")) %>%  # Make column of plot ID vectors 
  mutate(plot_id_length = sapply(.$plot_id_vec, length))  # Get number of plots in cluster

# Create plot ID / plot Cluster lookup table
plot_id_lookup <- plots %>% 
  filter(plot_id %in% unlist(plots_fil_sf$plot_id_vec)) %>%
  dplyr::select(plot_cluster, plot_id)

# Filter stems by plot ID
# Filter stems to >10 cm DBH
stem_size <- 10
stems_fil <- stems %>%
  inner_join(., plot_id_lookup, by = "plot_id") %>%
  filter(
    alive == "a",
    diam >= stem_size,
    plot_cluster %in% plots_fil_sf$plot_cluster)

# Summarise to only trees
trees <- stems_fil %>%
  filter(!is.na(species_name_clean)) %>%
  group_by(plot_id, tree_id) %>%
  summarise(
    plot_cluster = first(na.omit(plot_cluster)),
    species_name_clean = first(na.omit(species_name_clean)),
    diam = sum(diam, na.rm = TRUE),
    ba = sum(ba, na.rm = TRUE),
    agb = sum(agb, na.rm = TRUE))

# Find sites >5 species with >1 individual
plot_5sp <- trees %>% 
  group_by(plot_cluster, species_name_clean) %>%
  tally() %>%
  filter(n > 1) %>%
  group_by(plot_cluster) %>%
  tally() %>%
  filter(n > 5) %>%
  pull(plot_cluster)

# Find sites with 4 plots
plot_4p <- plot_id_lookup %>% 
  group_by(plot_cluster) %>% 
  tally() %>%
  filter(n == 4) %>%
  pull(plot_cluster)

# Find sites >50 trees ha
trees_ha <- 50
plot_50t <- trees %>% 
  filter(plot_cluster %in% plot_4p) %>%
  group_by(plot_cluster) %>%
  tally() %>%
  mutate(t_ha = n / 0.4) %>%
  filter(t_ha > trees_ha) %>%
  pull(plot_cluster)

# Find sites with no non-natives
plot_native <- trees %>% 
  group_by(plot_cluster) %>%
  summarise(exotics = ifelse(any(grepl("Pinus|Eucalyptus", species_name_clean)), TRUE, FALSE)) %>%
  filter(exotics == FALSE) %>% 
  pull(plot_cluster)

# Filter sites 
tree_fil <- trees %>%
  filter(plot_cluster %in% plot_5sp, 
    plot_cluster %in% plot_4p, 
    plot_cluster %in% plot_50t, 
    plot_cluster %in% plot_native)

# How many trees could not be identified?
ntrees <- length(tree_fil$species_name_clean)
nspindet <- sum(grepl("indet", tree_fil$species_name_clean) & 
  !grepl("Indet indet", tree_fil$species_name_clean)) 
perspindet <- nspindet / ntrees * 100
nuniquespindet <- length(unique(tree_fil$species_name_clean[grepl("indet", tree_fil$species_name_clean)]))
ngenindet <- sum(grepl("Indet indet", tree_fil$species_name_clean))
pergenindet <- ngenindet / ntrees * 100

# Find quadratic mean of tree diameter per site, and diam CoV
diam_summ <- tree_fil %>%
  group_by(plot_cluster) %>%
  summarise(
    diam_quad_mean = quadMean(diam, na.rm = TRUE),
    diam_mean = mean(diam, na.rm = TRUE),
    diam_sd = sd(diam, na.rm = TRUE)) %>%
  mutate(diam_cov = diam_sd / diam_mean * 100) %>% 
  dplyr::select(plot_cluster, diam_quad_mean, diam_cov)

# Create tree species abundance matrix by plot cluster
ba_clust_mat <- abMat(tree_fil, site_id = "plot_cluster", 
  species_id = "species_name_clean", abundance = "ba") %>%
  dplyr::select(-`Indet indet`)

# Filter plots data to match sites in filtered tree data
plots_clean <- plots_fil_sf %>% 
  filter(plot_cluster %in% tree_fil$plot_cluster) %>%
  left_join(., diam_summ, by = "plot_cluster")

# Are all plots in both objects?
stopifnot(nrow(ba_clust_mat) == nrow(plots_clean))

# Write files
saveRDS(plots_clean, "dat/plots.rds")
saveRDS(plot_id_lookup, "dat/plot_id_lookup.rds")
saveRDS(ba_clust_mat, "dat/ba_clust_mat.rds")

# Write stats to .tex
write(
  c(
    commandOutput(stem_size, "stemSize"),
    commandOutput(census, "censusDate"),
    commandOutput(n_total_sites, "nTotalSites"),
    commandOutput(ntrees, "nTrees"),
    commandOutput(trees_ha, "treesHa"),
    commandOutput(round(perspindet, 1), "perSpIndet"),
    commandOutput(round(pergenindet, 1), "perGenIndet")
    ),
  file = "out/data_prep_vars.tex")

