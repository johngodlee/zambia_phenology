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
plots <- read.csv("dat/plots_v2.7.csv")

stems <- read.csv("dat/stems_latest_v2.7.csv")

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

# Filter stems by plot ID
stem_size <- 10
stems_fil <- stems %>%
  inner_join(., plot_id_lookup, by = "plot_id") %>%
  filter(diam >= stem_size) %>% 
  filter(plot_cluster %in% plots_fil_sf$plot_cluster)

# Summarise to only trees
# Remove duplicated tree measurements
tree_fil <- stems_fil %>%
  filter(!is.na(species_name_clean)) %>%
  group_by(plot_id, tree_id) %>%
  summarise(
    plot_cluster = first(na.omit(plot_cluster)),
    species_name_clean = first(na.omit(species_name_clean)),
    diam = sum(diam, na.rm = TRUE),
    ba = sum(ba, na.rm = TRUE),
    agb = sum(agb, na.rm = TRUE))

# Create tree species abundance matrix by plot cluster
ba_clust_mat <- abMat(tree_fil, site_id = "plot_cluster", 
  species_id = "species_name_clean", abundance = "ba") 

# Occurrence based abundance matrix for filtering
occ_clust_mat <- abMat(tree_fil, site_id = "plot_cluster", 
  species_id = "species_name_clean", abundance = NULL)

# Remove genus level indets 
ba_clust_mat <- ba_clust_mat[,-which(names(ba_clust_mat) == "Indet indet")]

# Remove plots below 95th percentile of basal area 
ba_ha_lim <- quantile(plots_fil_sf$ba_ha, 0.05, na.rm = TRUE)

clust_plot_area <- plots_fil_sf$plot_id_length[match(rownames(ba_clust_mat), 
  plots_fil_sf$plot_cluster)] * 0.1
ba_clust_ba_ha <- rowSums(ba_clust_mat) / clust_plot_area
ba_ha_cluster <- row.names(ba_clust_mat[(ba_clust_ba_ha > ba_ha_lim),])
ba_clust_mat <- ba_clust_mat[row.names(ba_clust_mat) %in% ba_ha_cluster,]

# Remove plots with fewer than 5 species with more than 1 individual
sp_ab <- unname(apply(occ_clust_mat, 1, function(x) { 
    sum(x > 1, na.rm = TRUE) 
  }))
sp_ab_cluster <- row.names(occ_clust_mat[sp_ab >= 5,])
ba_clust_mat <- ba_clust_mat[row.names(ba_clust_mat) %in% sp_ab_cluster,] 

# Remove species with fewer than 5 total occurrences
occ_species <- names(occ_clust_mat[,colSums(occ_clust_mat) >= 5])
ba_clust_mat <- ba_clust_mat[,names(ba_clust_mat) %in% occ_species]

# Exclude species which are clearly non-native
ba_clust_mat <- ba_clust_mat[,-which(grepl("Pinus|Eucalyptus", names(ba_clust_mat)))] 

# Remove plots with no individuals
ba_clust_mat <- ba_clust_mat[rowSums(ba_clust_mat) != 0,]

# Clean other abundance matrices based on filtering above 
plots_clean <- plots_fil_sf[plots_fil_sf$plot_cluster %in% rownames(ba_clust_mat),]
occ_clust_mat_clean <- occ_clust_mat[row.names(occ_clust_mat) %in% row.names(ba_clust_mat), 
  names(occ_clust_mat) %in% names(ba_clust_mat)]

# Write files
saveRDS(plots_clean, "dat/plots.rds")
saveRDS(plot_id_lookup, "dat/plot_id_lookup.rds")
saveRDS(ba_clust_mat, "dat/ba_clust_mat.rds")
saveRDS(occ_clust_mat_clean, "dat/occ_clust_mat.rds")

# Write stats to .tex
write(
  c(
    commandOutput(ba_ha_lim, "baLim"),
    commandOutput(stem_size, "stemSize"),
    commandOutput(census, "censusDate"),
    commandOutput(n_total_sites, "nTotalSites")),
  file = "out/data_prep_vars.tex")
