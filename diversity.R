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
stems <- read.csv("dat/stems_latest_v2.7.csv")

plots <- readRDS("dat/plots_trmm.rds") 
plots <- st_as_sf(plots)

plot_id_lookup <- readRDS("dat/plot_id_lookup.rds")

# Filter stems by plot ID
stems_fil <- stems %>%
  inner_join(., plot_id_lookup, by = "plot_id") %>%
  filter(diam >= 10) %>% 
  filter(plot_cluster %in% plots$plot_cluster)

# Define stem percentage filter parameters
mopane_per <- 0.5
stems_ha <- 50
stem_size <- 10

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
      plots[plots$plot_cluster %in% row.names(.), "plot_id_length"]
      )) * 0.1) > stems_ha) %>%  # Remove plots with fewer than x stems ha
  filter((.$`Colophospermum mopane` / rowSums(.)) < mopane_per)  # Filter mopane plots

# Remove plots not in tree abundance matrix
plots_clean <- filter(plots, plot_cluster %in% rownames(tree_ab_mat_clust))
tree_ab_mat_plot_clean <- tree_ab_mat_plot[row.names(tree_ab_mat_plot) %in% 
  unlist(plots_clean$plot_id_vec),]

# Write tree abundance matrices
saveRDS(tree_ab_mat_clust, "dat/tree_ab_mat.rds")
saveRDS(tree_ab_mat_plot, "dat/tree_ab_mat_plot.rds")

# Exclude species which are clearly non-native
ab_mat_clean <- tree_ab_mat_clust %>%
  dplyr::select(-starts_with("Pinus"), -starts_with("Eucalyptus")) 

# Calculate common  diversity statistics
div_df <- data.frame(plot_cluster = row.names(ab_mat_clean), 
  richness = unname(rowSums(ab_mat_clean != 0)),
  shannon = diversity(ab_mat_clean),
  simpson = diversity(ab_mat_clean, "simpson"),
  evenness = diversity(ab_mat_clean) / log(rowSums(ab_mat_clean > 0)))

# Filter abundance matrix
# Remove plots with fewer than 5 species with more than 1 individual
ab_mat_fil <- ab_mat_clean[
  unname(apply(ab_mat_clean, 1, function(x) { sum(x > 1, na.rm = TRUE) })) >= 5,]  
ab_mat_fil <- ab_mat_fil[rowSums(ab_mat_fil) != 0,]  # Remove plots with no individuals
ab_mat_fil <- ab_mat_fil[,colSums(ab_mat_fil) > 5]  # Remove species with less than 5 occurrences

plots_fil <- plots_clean %>%
  filter(plot_cluster %in% row.names(ab_mat_fil))

# NSCA on species abundance
naxes <- 2
nsca <- dudi.nsc(df = ab_mat_fil, scannf = FALSE, nf = naxes)

nsca_inertia <- round(nsca$eig[naxes+1], 2)*10

# Extract euclidean distances between plots from NSCA
nsca_dist <- dist(nsca$li)

# Determine optimal number of clusters for hieriarchical clustering
kval <- seq(2,10)

v <- unlist(lapply(seq(length(kval)), function(x) {
    clust <- hclust(nsca_dist, method = "ward.D2")
    clust_cut <- cutree(clust, k = kval[x])
    ss <- cluster::silhouette(clust_cut, nsca_dist)
    mean(ss[, 3])
}))

sil <- data.frame(kval, v)

n_clusters <- sil$kval[which.max(sil$v)]

# Silhouette plot
pdf(file = "img/clust_sil.pdf", width = 8, height = 5)
ggplot(sil, aes(x = kval, y = v)) + 
  geom_point() + 
  geom_path() + 
  geom_vline(xintercept = n_clusters, linetype = 2) +
  theme_panel() + 
  labs(x = "Clusters", y = "Mean silhouette width") 
dev.off()

# Cluster euclidean distances of NSCA with ward algorithm
ward_clust <- hclust(nsca_dist, method = "ward.D2")
ward_dat <- dendro_data(ward_clust)

plot.new()
clust_classif <- rect.hclust(ward_clust, k = n_clusters)
clust_df <- do.call(rbind, lapply(seq(length(clust_classif)), function(x) {
  data.frame(plot_cluster = names(clust_classif[[x]]), cluster = as.character(x))
}))

ward_merge <- left_join(ward_dat$labels, clust_df, by = c("label" = "plot_cluster"))

ward_rect <- ward_merge %>% 
  group_by(cluster) %>%
  summarise(xmin = min(x)+1, 
    xmax = max(x)-1) %>%
  mutate(ymin = 0, ymax = 3,
    cluster = factor(cluster, 
      labels = clust_lookup[1:length(unique(.$cluster))]))

# Dendrogram with clusters
pdf(file = "img/clust_dendro.pdf", width = 8, height = 5)
ggplot() + 
  geom_segment(data = ward_dat$segments, 
    aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_rect(data = ward_rect, 
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = cluster),
    size = 1.2, fill = NA) + 
  theme_panel() + 
  scale_colour_manual(values = clust_pal)
dev.off()

# Extract data from clustering and NSCA
nsc_df <- ward_merge %>%
  dplyr::select(plot_cluster = label, cluster) %>%
  left_join(., rownames_to_column(nsca$l1), by = c("plot_cluster" = "rowname")) %>%
  mutate(cluster = as.character(cluster))

# NSCA ordination plot with clusters
nscaPlot <- function(x,y, dat = nsc_df) {
  x <- ensym(x)
  y <- ensym(y)

  p <- ggplot(dat, aes(x = !!x, y = !!y)) +
    geom_point(shape = 21, size = 2, aes(fill = cluster)) +
    geom_bag(prop = 0.95, alpha = 0.2, aes(colour = cluster, fill = cluster)) +
    scale_fill_manual(name = "Cluster", values = clust_pal) +
    scale_colour_manual(name = "Cluster", values = clust_pal) +
    theme_panel() +
    labs(x = gsub("RS", "NSC ", x),
      y = gsub("RS", "NSC ", y))

  return(p)
}

pdf(file = "img/nsca.pdf", width = 8, height = 6)
nscaPlot(RS1, RS2) +
  plot_layout(guides = "collect", widths = 1) &
  scale_colour_manual(name = "Cluster", values = brightness(clust_pal, 0.5),
    limits = unique(nsc_df$cluster))
dev.off()

# Add values to data
div <- plots_fil %>%
  left_join(., div_df, by = "plot_cluster")  %>%
  left_join(., clust_df, by = "plot_cluster") 

# Write file
saveRDS(div, "dat/plots_div.rds") 

# Indicator species per cluster
clust_indval <- indval(ab_mat_fil, clustering = div$cluster)

# Summarise indicator analysis
summary(clust_indval, p = 0.05, type = "short", digits = 2, show = p)
clust_indval$indval$sp <- row.names(clust_indval$indval)
row.names(clust_indval$indval) <- seq(from = 1, to = length(clust_indval$indval$sp))

indval_extrac <- lapply(1:n_clusters, function(x) {
    out <- head(clust_indval$indval[order(clust_indval$indval[[x]], 
          decreasing = TRUE),c(n_clusters+1, x)])
    out[!grepl("indet", out$sp, ignore.case = TRUE),]
  })

indval_extrac_tidy <- do.call(rbind, lapply(indval_extrac, function(x) {
    cluster <- names(x)[2]
    out <- x[1:3,]
    out$cluster <- cluster
    names(out) <- c("species", "indval", "cluster")
    out[,c(3,1,2)]
  })
)

# Get species richness for each cluster +/- IQR
# Get number of sites
# Get MAP mean and Diurnal dT mean
clust_summ <- div %>%
  st_drop_geometry()%>%
  group_by(cluster) %>%
  summarise(richness_median = median(richness, na.rm = TRUE),
    richness_iqr = (quantile(richness, 0.75) - quantile(richness, 0.25)),
    n_sites = as.character(n()),
    map_mean = mean(map, na.rm = TRUE),
    map_sd = sd(map, na.rm = TRUE),
    diurnal_dt_mean = mean(diurnal_temp_range, na.rm = TRUE),
    diurnal_dt_sd = sd(diurnal_temp_range, na.rm = TRUE)) %>%
  mutate(richness = paste0(sprintf("%.0f", richness_median), "(", sprintf("%.0f", richness_iqr), ")"),
    map = paste0(sprintf("%.0f", map_mean), "(", sprintf("%.1f", map_sd), ")"),
    diurnal_dt = paste0(sprintf("%.0f", diurnal_dt_mean), "(", sprintf("%.1f", diurnal_dt_sd), ")")) %>%
  dplyr::select(cluster, n_sites, richness, map, diurnal_dt) %>%
  left_join(., indval_extrac_tidy, by = "cluster") %>%
  mutate(species = paste0("\\textit{", species, "}"),
  indval = sprintf("%.3f", indval))

clust_summ[c(1,3,4,6,7,9),c(1,2,3,4,5)] <- ""

names(clust_summ) <- c("Cluster", "N sites", "Richness", "MAP", "Diurnal $\\delta$T", "Species", "Indicator value")

# Export indval table
clust_summ_xtable <- xtable(clust_summ, 
  label = "clust_summ",
  align = rep("c", 8),
  display = rep("s", 8),
  caption = "Climatic information and Dufrene-Legendre indicator species analysis for the vegetation type clusters identified by the PAM algorithm. The three species per cluster with the highest indicator values are shown along with other key statistics for each cluster. MAP (Mean Annual Precipitation) and Diurnal $\\delta$T are reported as the mean and 1 standard deviation in parentheses. Species richness is reported as the median and the interquartile range in parentheses.")

fileConn <- file("out/clust_summ.tex")
writeLines(print(clust_summ_xtable, include.rownames = FALSE,
    table.placement = "H",
    hline.after = c(-1,0,seq(from = 3, by = 3, length.out = n_clusters-1)),
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Test Beta diversity between plots in a cluster - UNFINISHED

# Filter to plots we're using
tree_ab_mat_plot_fil <- tree_ab_mat_plot[row.names(tree_ab_mat_plot) %in% 
  unlist(strsplit(div$plot_id, ",")) & rowSums(tree_ab_mat_plot) != 0,]

# Split each plot cluster into own abundance matrix
tree_ab_mat_split <- left_join(rownames_to_column(tree_ab_mat_plot_fil), plot_id_lookup, 
  by = c("rowname" = "plot_id")) %>%
  dplyr::select(-rowname) %>%
  split(., .$plot_cluster)

# Clean each abundance matrix to remove columns with 0 or 1 record
tree_ab_split_clean <- lapply(tree_ab_mat_split, function(x) {
  x %>%
    dplyr::select(-plot_cluster) %>%
    dplyr::select(which(!colSums(., na.rm = TRUE) %in% c(0, 1)))
  })

# Calculate mean pairwise Bray distances within each cluster
plot_dist_mean <- unlist(lapply(tree_ab_split_clean, function(x) { 
  mean(vegdist(x))
  }))

# Calculate mean pairwise Bray distance across all plots
plot_dist_all_mean <- mean(vegdist(tree_ab_mat_plot_fil))

# Create histogram
pdf(file = "img/plot_dist_hist.pdf", width = 8, height = 6)
ggplot() + 
  geom_histogram(data = data.frame(plot_dist_mean), aes(x = plot_dist_mean), 
    fill = pal[6], colour = "black", binwidth = 0.05 ) +
  geom_vline(xintercept = plot_dist_all_mean, colour = "red") + 
  theme_panel() + 
  labs(x = "Mean pairwise Bray distance within cluster", y = "Frequency")
dev.off()

# What percentage of clusters have a mean pairwise distance lower than the mean across all pairs? 
plot_dist_mean_clean <- plot_dist_mean[!is.na(plot_dist_mean)]
plot_dist_per <- round(length(which(plot_dist_mean_clean < plot_dist_all_mean)) / 
  length(plot_dist_mean_clean) * 100, 1)

# How many plots are in each cluster?
clust_tally <- div %>% 
  count(cluster) %>% 
  st_drop_geometry() %>%
  pull(n)

write(
  c(
    commandOutput(plot_dist_per, "plotDistPer"),
    commandOutput(mopane_per*100, "mopanePer"),
    commandOutput(nsca_inertia, "nscaInertia"),
    commandOutput(naxes, "nscaAxes"),
    commandOutput(stems_ha, "stemsHa"),
    commandOutput(stem_size, "stemSize"),
    commandOutput(n_clusters, "nCluster"),
    commandOutput(clust_tally[1], "nClusterA"),
    commandOutput(clust_tally[2], "nClusterB"),
    commandOutput(clust_tally[3], "nClusterC")
    ),
  file = "out/diversity_vars.tex")

