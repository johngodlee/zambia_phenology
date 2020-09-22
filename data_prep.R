# Filter SEOSAW data
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(dplyr)
library(tidyr)
library(sf)

source("functions.R")

# Delete variables TeX file 
if (file.exists("out/vars.tex")) {
  file.remove("out/vars.tex")
}

# Import data
plots <- read.csv("~/git_proj/seosaw_data/data_out/plots_v2.7.csv")
stems <- read.csv("~/git_proj/seosaw_data/data_out/stems_latest_v2.7.csv")
load("/Volumes/john/seosaw_plot_summary5Apr2019.Rdata")


# Create clean clusters dataframe
clust_df <- ssaw8$cluster %>%
  filter(file == "46_zambia_siampale_ILUAii_checkN.xlsx") %>%
  dplyr::select(plot_id = plotcode, contains("clust"))

# Subset plots and data fields - Zambia ILUAii data only
plots_fil_sf <- plots %>% 
  filter(
    prinv == "Siampale A.") %>% 
  dplyr::select(
    plot_id, plot_cluster, longitude_of_centre, latitude_of_centre, 
    richness, shannon, simpson, evenness, n_stems_gt5_ha, ba_ha, agb_ha,
    mat = bio1, diurnal_temp_range = bio2, map = bio12,
    clay = CLYPPT, sand = SNDPPT, cec = CECSOL,
    last_census_date) %>%
  left_join(., clust_df, by = "plot_id") %>%
  group_by(plot_cluster) %>%
  summarise(
    plot_id = paste0(plot_id, collapse = ","),
    clust7 = first(na.omit(clust7)),
    clust5 = first(na.omit(clust5)),
    clust4 = first(na.omit(clust4)),
    longitude_of_centre = mean(longitude_of_centre, na.rm = TRUE),
    latitude_of_centre = mean(latitude_of_centre, na.rm = TRUE),
    n_stems_gt5_ha = mean(n_stems_gt5_ha, na.rm = TRUE),
    ba_ha = mean(ba_ha, na.rm = TRUE),
    agb_ha = mean(agb_ha, na.rm = TRUE),
    mat = mean(mat, na.rm = TRUE),
    diurnal_temp_range = mean(diurnal_temp_range, na.rm = TRUE),
    map = mean(map, na.rm = TRUE),
    clay = mean(clay, na.rm = TRUE),
    sand = mean(sand, na.rm = TRUE),
    cec = mean(cec, na.rm = TRUE),
    last_census_date = first(last_census_date)) %>%
  st_as_sf(., coords = c("longitude_of_centre", "latitude_of_centre")) %>%
  `st_crs<-`(4326) %>%
  mutate(plot_id_vec = strsplit(as.character(plot_id), 
  split = ",")) %>%
  mutate(plot_id_length = sapply(.$plot_id_vec, length)) %>%
  filter(plot_id_length == 4,
    !is.na(clust5)) %>%
  dplyr::select(-plot_id_length)

census <- unique(plots_fil_sf$last_census_date)

n_total_sites <- plots %>%
  filter(prinv == "Siampale A.") %>%
  pull(plot_cluster) %>%
  unique() %>%
  length()

write(
  c(
    commandOutput(census, "censusDate"),
    commandOutput(n_total_sites, "nTotalSites")
    ),
  file="out/vars.tex", append=TRUE)


plots_fil_sf <- plots_fil_sf %>% 
  dplyr::select(-last_census_date)

# Create plot ID / plot Cluster lookup table
plot_id_lookup <- plots %>% 
  filter(plot_id %in% unlist(plots_fil_sf$plot_id_vec)) %>%
  dplyr::select(plot_cluster, plot_id)

saveRDS(plot_id_lookup, "dat/plot_id_lookup.rds")

# Filter stems by plots
stems_fil <- stems %>%
  inner_join(., plot_id_lookup, by = "plot_id") %>%
  dplyr::select(-plot_id)

# Mopane stem percentage filter
mopane_per <- 0.5
stems_ha <- 50

write(
  c(
    commandOutput(stems_ha, "stemsHa"),
    commandOutput(mopane_per, "mopanePer")
  ),
  file="out/vars.tex", append=TRUE)

# Create tree species abundance matrix
tree_ab_mat <- stems_fil %>% 
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
  dplyr::select(-`Indet indet`)

# Clean and filter tree species abundance matrix
rownames(tree_ab_mat) <- tree_ab_mat[["plot_cluster"]]
tree_ab_mat <- tree_ab_mat[,-1, drop=FALSE]  # Remove plot cluster column
tree_ab_mat <- filter(tree_ab_mat, (tree_ab_mat$`Colophospermum mopane` / rowSums(tree_ab_mat)) < mopane_per)  # Filter mopane plots
tree_ab_mat <- filter_all(tree_ab_mat, any_vars(. != 0))  # Filter empty plots
tree_ab_mat <- tree_ab_mat[stems_ha < rowSums(tree_ab_mat) / 0.4,]  # Remove plots with fewer than x stems ha

# Remove plots not in tree abundance matrix
plots_fil_sf <- filter(plots_fil_sf, plot_cluster %in% rownames(tree_ab_mat))

saveRDS(plots_fil_sf, "dat/plots.rds")
saveRDS(tree_ab_mat, "dat/tree_ab_mat.rds")

