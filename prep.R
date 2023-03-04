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
library(sf)

source("tex_func.R")

# Import data
plots <- read.csv("dat/plots_v2.12.csv")
stems <- read.csv("dat/stems_iluaii_v2.12.csv")

# Subset plots to Zambia ILUAii data only
plots_iluaii <- plots %>% 
  filter(
    grepl("ZIS", plot_id), 
    prinv == "Siampale A.")

# Find census year of the data
census <- unique(plots_iluaii$census_date) 

# Find number of sites (not plots)
n_total_sites <- plots_iluaii %>%
  pull(plot_cluster) %>%
  unique() %>%
  length()

# Summarise plots to sites and select columns
plots_fil_sf <- plots_iluaii %>%  
  dplyr::select(
    plot_id,
    plot_cluster,
    longitude_of_centre,
    latitude_of_centre) %>%  # Select columns
  group_by(plot_cluster) %>%  # Group by plot cluster
  summarise(
    plot_id = paste0(plot_id, collapse = ","),
    longitude_of_centre = mean(longitude_of_centre, na.rm = TRUE),
    latitude_of_centre = mean(latitude_of_centre, na.rm = TRUE)) %>%  # Summarise by plot cluster
  st_as_sf(., coords = c("longitude_of_centre", "latitude_of_centre"), 
    crs = 4326) %>%  # Make spatial
  mutate(plot_id_vec = strsplit(as.character(plot_id), split = ","))  # Plot ID vector
    
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

# Summarise stems to trees
trees <- stems_fil %>%
  filter(!is.na(species_name_clean)) %>%
  group_by(plot_id, tree_id) %>%
  summarise(
    plot_cluster = first(na.omit(plot_cluster)),
    species_name_clean = first(na.omit(species_name_clean)),
    diam = sum(diam, na.rm = TRUE),
    ba = sum(ba, na.rm = TRUE))

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

# Find sites with trees in at least 2 plots
site_2wt <- trees %>% 
  group_by(plot_cluster, plot_id) %>% 
  tally() %>% 
  mutate(n = ifelse(n >0, 1, 0)) %>% 
  group_by(plot_cluster) %>% 
  tally() %>% 
  filter(n > 1) %>% 
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
  summarise(exotics = ifelse(
      any(grepl("Pinus|Eucalyptus", species_name_clean)), TRUE, FALSE)) %>%
  filter(exotics == FALSE) %>% 
  pull(plot_cluster)

# Exclude sites based on filters above
tree_fil <- trees %>%
  filter(
    plot_cluster %in% plot_5sp, 
    plot_cluster %in% plot_4p, 
    plot_cluster %in% plot_50t, 
    plot_cluster %in% site_2wt,
    plot_cluster %in% plot_native)

# Filter plots data to match sites in filtered tree data
plots_clean <- plots_fil_sf %>% 
  filter(plot_cluster %in% tree_fil$plot_cluster)

# Check all sites in both dataframes
stopifnot(all(sort(unique(tree_fil$plot_cluster)) == sort(plots_clean$plot_cluster)))

# Write files
saveRDS(plots_clean, "dat/plots.rds")
saveRDS(tree_fil, "dat/trees.rds")

# Write stats to .tex
write(
  c(
    commandOutput(stem_size, "stemSize"),
    commandOutput(census, "censusDate"),
    commandOutput(n_total_sites, "nTotalSites"),
    commandOutput(trees_ha, "treesHa")

    ),
  file = "out/prep_vars.tex")

