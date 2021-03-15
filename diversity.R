# Get diversity statistics
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(dplyr)
library(tidyr)
library(sf)
library(vegan)
library(ade4)
library(ggplot2)
library(patchwork)
library(cluster)
library(labdsv)
library(shades)
library(tibble)
library(xtable)
library(ggdendro)

source("functions.R")

# Import data
plots <- readRDS("dat/plots_try.rds") 
plots <- st_as_sf(plots)

plot_id_lookup <- readRDS("dat/plot_id_lookup.rds")

ba_clust_mat <- readRDS("dat/ba_clust_mat.rds")

# Calculate common diversity statistics
div_df <- data.frame(plot_cluster = row.names(ba_clust_mat), 
  richness = unname(rowSums(ba_clust_mat != 0)),
  shannon = diversity(ba_clust_mat),
  simpson = diversity(ba_clust_mat, "simpson"),
  evenness = diversity(ba_clust_mat) / log(rowSums(ba_clust_mat > 0)))

# Effective true numbers diversity
# based on shannon calculated by weighting on basal area.
div_df$eff_rich <- exp(div_df$shannon)

# Gather matrix to dataframe
ba_gather <- ba_clust_mat %>%
  rownames_to_column("plot_cluster") %>%
  gather(species, ba, -plot_cluster) %>%
  filter(ba > 0) %>%
  group_by(plot_cluster, species) %>%
  summarise(ba = sum(ba, na.rm = TRUE)) 

# Split by plot_cluster
ba_split <- split(ba_gather, ba_gather$plot_cluster)

# Dominant species per plot by basal area
dom_sp <- do.call(rbind, lapply(ba_split, function(x) {
  ba_total <- sum(x$ba, na.rm = TRUE)
  x_order <- x[order(x$ba, decreasing = TRUE),]
  x_fil <- x_order[1:5,]
  x_fil$prop <- x_fil$ba / ba_total
  x_fil$id <- 1:5
  x_spread <- pivot_wider(x_fil, names_from = id, 
    values_from = c("species", "ba", "prop"))

  return(x_spread)
}))

# Add values to data
# How many of top n dominant species make up percentages of basal area?
plots_div <- plots %>%
  left_join(., div_df, by = "plot_cluster") %>%
  left_join(., dom_sp, by = "plot_cluster") %>%
  mutate(
    prop_1_cum = prop_1,
    prop_2_cum = prop_1_cum + prop_2,
    prop_3_cum = prop_2_cum + prop_3,
    prop_4_cum = prop_3_cum + prop_4,
    prop_5_cum = prop_4_cum + prop_5)

plots_div %>%
  st_drop_geometry() %>%
  dplyr::select(ends_with("_cum")) %>%
  gather(key, value) %>%
  ggplot(., aes(x = value)) + 
  geom_histogram(colour = "black", fill = "grey") + 
  geom_vline(xintercept = 0.8, linetype = 2, colour = "red") + 
  facet_wrap(~key)

# Reverse some columns so all increase with conservativism
rev_cols <- c("leaf_n_mass_genus_cwm", "leaf_p_mass_genus_cwm", "sla_genus_cwm")
cons_acq <- plots_div %>%
  dplyr::select(
    plot_cluster, 
    ends_with("_genus_cwm"), 
    -wood_n_mass_genus_cwm,
    -leaf_thick_genus_cwm, 
    -bark_thick_genus_cwm,
    -leaf_n_area_genus_cwm) %>%
  mutate(across(all_of(rev_cols), 
      ~-1 * .x,
      .names = "{.col}_rev")) %>%
  dplyr::select(-all_of(rev_cols)) %>%
  filter(across(everything(), ~!is.na(.x)))

# PCA of community weighted means
cons_pca <- prcomp(st_drop_geometry(cons_acq[,-which(names(cons_acq) == "plot_cluster")]), 
  center = TRUE, scale. = TRUE)

cons_pca_tidy <- as.data.frame(cons_pca$x) %>%
  mutate(plot_cluster = cons_acq$plot_cluster)

cons_pca_df <- plots_div %>%
  left_join(., cons_pca_tidy, by = "plot_cluster") %>%
  mutate(genus_1 = unlist(lapply(strsplit(species_1, " "), "[", 1)))

ggplot() + 
  geom_point(data = cons_pca_df, 
    aes(x = PC1, y = PC2, fill = genus_1, size = prop_1),
    shape = 21, colour = "black")

# Write file
saveRDS(plots_div, "dat/plots_div.rds") 
