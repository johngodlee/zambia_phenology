# Filter SEOSAW data for phenology analysis
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(dplyr)
library(tidyr)
library(sf)
library(tibble)

source("functions.R")

# Delete variables TeX file 
if (file.exists("out/vars.tex")) {
  file.remove("out/vars.tex")
}

# Import data
plots <- read.csv("~/git_proj/seosaw_data/data_out/plots_v2.7.csv")
stems <- read.csv("~/git_proj/seosaw_data/data_out/stems_latest_v2.7.csv")
load("/Users/johngodlee/google_drive/phd/thesis/regional_befr/data/seosaw_plot_summary5Apr2019.Rdata")

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
stems_fil <- stems %>%
  inner_join(., plot_id_lookup, by = "plot_id") %>%
  filter(diam >= 10)

# Define stem percentage filter parameters
mopane_per <- 0.5
stems_ha <- 50
stem_size <- 10

# Create tree species abundance matrix by plot
tree_ab_mat_plot <- stems_fil %>% 
  dplyr::select(plot_id, tree_id, species_name_clean) %>%  # Select columns
  filter(!is.na(species_name_clean)) %>%  # Remove stems with no species
  group_by(plot_id, tree_id) %>%  # Group by plot and tree ID
  filter(row_number() == 1) %>%  # Remove duplicated tree measurements
  group_by(plot_id, species_name_clean, .drop = FALSE) %>%
  tally() %>%
  spread(species_name_clean, n, fill = 0) %>%
  ungroup() %>%
  mutate_at(vars(-plot_id), as.double) %>%
  as.data.frame() %>%
  dplyr::select(-`Indet indet`) %>%  # Remove stems with no species
  column_to_rownames("plot_id")

# Create tree species abundance matrix by plot cluster
tree_ab_mat_clust <- stems_fil %>% 
  dplyr::select(plot_cluster, tree_id, species_name_clean) %>%
  filter(!is.na(species_name_clean)) %>%
  group_by(plot_cluster, tree_id) %>%
  filter(row_number() == 1) %>%
  group_by(plot_cluster, species_name_clean, .drop = FALSE) %>%
  tally() %>%
  spread(species_name_clean, n, fill = 0) %>%
  ungroup() %>%
  mutate_at(vars(-plot_cluster), as.double) %>%
  as.data.frame() %>%
  dplyr::select(-`Indet indet`) %>%
  column_to_rownames("plot_cluster") %>%
  filter_all(any_vars(. != 0)) %>%
  filter(rowSums(.) / 
    (pull(st_drop_geometry(
      plots_fil_sf[plots_fil_sf$plot_cluster %in% row.names(.), "plot_id_length"]
      )) * 0.1) > stems_ha) %>%  # Remove plots with fewer than x stems ha
  filter((.$`Colophospermum mopane` / rowSums(.)) < mopane_per)  # Filter mopane plots

# Remove plots not in tree abundance matrix
plots_clean <- filter(plots_fil_sf, plot_cluster %in% rownames(tree_ab_mat_clust))

tree_ab_mat_plot_clean <- tree_ab_mat_plot[row.names(tree_ab_mat_plot) %in% 
  unlist(plots_clean$plot_id_vec),]

# Write files
saveRDS(plots_clean, "dat/plots.rds")
saveRDS(tree_ab_mat_clust, "dat/tree_ab_mat.rds")
saveRDS(tree_ab_mat_plot, "dat/tree_ab_mat_plot.rds")
saveRDS(plot_id_lookup, "dat/plot_id_lookup.rds")

# Write stats to .tex
write(
  c(commandOutput(census, "censusDate"),
    commandOutput(stems_ha, "stemsHa"),
    commandOutput(stem_size, "stemSize"),
    commandOutput(mopane_per*100, "mopanePer"),
    commandOutput(n_total_sites, "nTotalSites")),
  file="out/vars.tex", append=TRUE)
