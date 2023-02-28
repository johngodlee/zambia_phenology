# Get diversity statistics
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(sf)
library(vegan)
library(cluster)
library(labdsv)
library(taxize)

source("plot_func.R")
source("tex_func.R")

# Import data
plots <- readRDS("dat/plots.rds") 
trees <- readRDS("dat/trees.rds")
ab_plot_mat <- readRDS("dat/ab_plot_mat.rds")
ba_clust_mat <- readRDS("dat/ba_clust_mat.rds")

# Get taxonomic ranks of all species
sp_tax <- gnr_resolve(names(ba_clust_mat),
  resolve_once = FALSE, with_context = TRUE, 
  canonical = TRUE, with_canonical_ranks = TRUE,
  best_match_only = TRUE, cap_first = FALSE, fields = "all",
  data_source_ids = c(179, 167, 165, 11, 4))

# Extract family and subfamily
split_rank <- strsplit(sp_tax$classification_path_ranks, split = "\\|")

fam_locs <- unlist(lapply(split_rank, function(x) {
  if (length(which(x == "family")) == 1) {
    which(x == "family")
  } else {
    NA
  }
}))

subfam_locs <- unlist(lapply(split_rank, function(x) {
  if (length(which(x == "subfamily")) == 1) {
    which(x == "subfamily")
  } else {
    NA
  }
}))

split_taxon <- strsplit(sp_tax$classification_path, split = "\\|")

sp_tax$family <- unlist(lapply(1:length(split_taxon), function(x) {
  gsub("^$", NA, split_taxon[[x]][fam_locs[x]])
}))

sp_tax$subfamily <- unlist(lapply(1:length(split_taxon), function(x) {
  gsub("^$", NA, split_taxon[[x]][subfam_locs[x]])
}))

# Write species information to file
saveRDS(sp_tax, "dat/taxon.rds")

# Calculate common diversity statistics
div <- data.frame(
  plot_cluster = row.names(ba_clust_mat), 
  richness = unname(rowSums(ba_clust_mat != 0)),
  shannon = diversity(ba_clust_mat),
  simpson = diversity(ba_clust_mat, "simpson"),
  evenness = diversity(ba_clust_mat) / log(rowSums(ba_clust_mat > 0)))

# Effective true numbers diversity
# based on Shannon calculated by weighting on basal area.
div$eff_rich <- exp(div$shannon)

# Check all plots included in diversity statistics
stopifnot(nrow(div) == nrow(plots))

# Gather basal area matrix to dataframe
ba_gather <- ba_clust_mat %>%
  rownames_to_column("plot_cluster") %>%
  gather(species, ba, -plot_cluster) %>%
  filter(ba > 0) %>%
  group_by(plot_cluster, species) %>%
  summarise(ba = sum(ba, na.rm = TRUE)) 

# Find percentage of ba which is Detarioideae and other families / subfamilies
ba_perc <- ba_gather %>% 
  left_join(., sp_tax[,c("user_supplied_name", "family", "subfamily")], 
    by = c("species" = "user_supplied_name")) %>%
  group_by(plot_cluster) %>%
  mutate(ba_total = sum(ba, na.rm = TRUE)) %>%
  mutate(ba_perc = ba / ba_total) 

ba_perc_fam <- ba_perc %>%
  group_by(plot_cluster, family) %>%
  summarise(ba_perc_fam = sum(ba_perc, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(plot_cluster, family, fill = list(ba_perc_fam = 0)) 

ba_perc_subfab_wide <- ba_perc %>%
  filter(family == "Fabaceae") %>% 
  filter(!is.na(subfamily)) %>%
  group_by(plot_cluster, subfamily) %>%
  summarise(ba_perc_subfab = sum(ba_perc, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(plot_cluster, subfamily, fill = list(ba_perc_subfab = 0)) %>%
  spread(subfamily, ba_perc_subfab) %>%
  mutate(across(everything(), ~ifelse(is.na(.x), 0, .x)) )

# Check no sites are duplicated
stopifnot(all(!duplicated(ba_perc_subfab_wide$plot_cluster)))

# For each plot_cluster, dominant species and genera by basal area
dom_sp <- ba_gather %>% 
  group_by(plot_cluster) %>%
  mutate(species_prop = ba / sum(ba, na.rm = TRUE)) %>%
  slice_max(order_by = ba, n = 10, with_ties = FALSE) %>%
  rename(species_ba = ba ) %>%
  mutate(id = row_number()) %>%
  left_join(., sp_tax[,c("user_supplied_name", "family", "subfamily")], 
    by = c("species" = "user_supplied_name")) %>%
  rename(species_family = family, species_subfamily = subfamily) %>%
  mutate(species_prop_cum = cumsum(species_prop)) %>%
  pivot_wider(names_from = id,
    values_from = c("species", "species_ba", "species_prop", 
      "species_family", "species_subfamily", "species_prop_cum")) %>%
  ungroup()

# Check all sites included
stopifnot(nrow(ba_clust_mat) == length(unique(dom_sp$plot_cluster)))
stopifnot(nrow(ba_clust_mat) == nrow(plots))

ba_chosen_mat <- ba_clust_mat

# Construct distance matrix
ba_dist <- vegdist(ba_chosen_mat, method = "bray")

# Ward distance
beta_ward <- hclust(ba_dist, method = "ward.D2")

# For each number of clusters create dataframe of dominant species
nclust <- seq(2, 25, 1)
dom_all <- lapply(nclust, function(y) {

  cuts <- cutree(beta_ward, y)

  sil <- silhouette(cuts, ba_dist)

  grp <- data.frame(
    plot_cluster = row.names(ba_chosen_mat),
    cu = cuts)
  names(grp)[2] <- paste0("c", y)

  tmp <- aggregate(ba_chosen_mat, list(grp[,2]), FUN = mean)
  dom <- t(round(tmp[,-1],1))

  dom_df <- do.call(rbind, lapply(seq(1, y), function(z) {
    ret <- sort(dom[,z], decreasing = TRUE)[1:5]
    data.frame(y, z, seq(1,5), names(ret), unname(ret))
    }))
  names(dom_df) <- c("nclust", "clust", "ndom", "taxa", "ab")

  list(grp, dom_df, sil)
})

# Extract silhouette 
sil_out <- lapply(dom_all, "[[", 3)

# Get mean silhouette value for each number of clusters
sil_mean <- unlist(lapply(sil_out, function(x) { 
  if (is.matrix(x)) {
    mean(x[,3], na.rm = TRUE)
  } else { 
    NA_real_
  }
}))

sil_mean_best <- max(sil_mean) * 10

# Identify optimal number of clusters
clust_optim <- nclust[sil_mean == max(sil_mean) & !is.na(sil_mean)]
##' 9
clust_optim <- 4

# Extract silhouette values for optimal cluster number
sil_best <- as.data.frame(matrix(c(sil_out[[clust_optim]]), ncol = 3))

names(sil_best) <- c("cluster", "neighbor", "sil_width")

sil_best_order <- sil_best[rev(order(sil_best$sil_width)),]

sil_best_order$id <- seq_len(nrow(sil_best_order))

# Plot silhouette values per cluster for optimal number of clusters
pdf(file = "img/ward_sil.pdf", width = 6, height = 4)
ggplot() + 
  geom_bar(data = sil_best_order, stat = "identity", position = "dodge",
    aes(x = cluster, y = sil_width, fill = as.character(cluster), group = id)) +
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  theme_bw() + 
  labs(x = "Cluster", y = "Silhouette width")
dev.off()

# Plot mean silhouette for each number of clusters
pdf(file = "img/ward_sil_mean.pdf", width = 6, height = 4)
ggplot() + 
  geom_point(aes(y = sil_mean, x = nclust)) + 
  geom_line(aes(y = sil_mean, x = nclust)) + 
  geom_vline(xintercept = clust_optim, linetype = 2) + 
  theme_bw() + 
  labs(x = "N clusters", y = "Mean silhouette width")
dev.off()

# Create dataframe of dominant species per cluster for numbers of clusters
dom_out <- do.call(rbind, lapply(dom_all, "[[", 2))

# Extract cluster of each plot per number of clusters
clust_out <- Reduce( function(...) { merge(..., by = "plot_cluster") }, 
  lapply(dom_all, "[[", 1) )

clust_best <- clust_out[,c("plot_cluster", paste0("c", clust_optim))]

# Indicator species per cluster
clust_indval <- indval(ba_chosen_mat, clustering = clust_best[[paste0("c", clust_optim)]])

# Summarise indicator analysis
clust_indval$indval$sp <- row.names(clust_indval$indval)
row.names(clust_indval$indval) <- seq(from = 1, to = length(clust_indval$indval$sp))

indval_extrac <- lapply(1:clust_optim, function(x) {
    out <- head(clust_indval$indval[
      order(clust_indval$indval[[x]], decreasing = TRUE),
      c(clust_optim+1, x)])
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

saveRDS(indval_extrac_tidy, "./dat/indval.rds")

# How many of top n dominant species make up percentages of basal area?
pdf(file = "img/basal_area_dom_hist.pdf", width = 10, height = 8)
dom_sp %>%
  dplyr::select(plot_cluster, contains("_cum_")) %>%
  gather(key, value, -plot_cluster) %>%
  filter(!is.na(value)) %>% 
  ggplot(., aes(x = value)) + 
  geom_histogram(colour = "black", fill = "grey") + 
  geom_vline(xintercept = 0.8, linetype = 2, colour = "red") + 
  facet_wrap(~key) + 
  theme_panel()
dev.off()

# Find quadratic mean of tree diameter per site
diam_summ <- trees %>%
  group_by(plot_cluster) %>%
  summarise(diam_quad_mean = sqrt(mean(diam^2, na.rm = TRUE))) %>%
  dplyr::select(plot_cluster, diam_quad_mean)

# Add values to data
plots_div <- plots %>%
  dplyr::select(plot_cluster) %>% 
  st_drop_geometry() %>% 
  left_join(., div, by = "plot_cluster") %>%  # Diversity stats
  left_join(., ba_perc_subfab_wide, by = "plot_cluster") %>%  # Dominant species
  left_join(., dom_sp, by = "plot_cluster") %>%  # Dominant species
  left_join(., clust_best, by = "plot_cluster") %>%  # Ward clusters 
  left_join(., diam_summ, by = "plot_cluster") %>%  # Quadratic mean
  rename(cluster = paste0("c", clust_optim))

# Check no duplications or removals
stopifnot(nrow(plots) == nrow(plots_div))
stopifnot(all(!is.na(plots_div$plot_cluster)))

# Calculate mean pairwise Bray distances within each cluster
# Filter to plots we're using 
plot_id_spread <- plots %>% 
  st_drop_geometry() %>% 
  dplyr::select(plot_cluster, plot_id) %>% 
  separate_rows(plot_id, sep = ",")

tree_ab_mat_plot_fil <- ab_plot_mat %>%
  rownames_to_column("plot_id") %>%
  left_join(., plot_id_spread, by = "plot_id") %>%
  dplyr::select(-plot_id) 

# Mean pairwise Bray distance across all plots
tree_ab_mat_plot_fil_fil <- tree_ab_mat_plot_fil %>%
  dplyr::select(-plot_cluster) %>%
  filter(rowSums(.) != 0)
plot_dist_all_mean <- mean(vegdist(tree_ab_mat_plot_fil_fil))

# Split by plot cluster
tree_ab_mat_plot_fil_split <- split(tree_ab_mat_plot_fil, tree_ab_mat_plot_fil$plot_cluster)

stopifnot(length(tree_ab_mat_plot_fil_split) == nrow(plots))

# Calculate Bray distance among plots in each cluster
plot_dist_mean <- unlist(lapply(tree_ab_mat_plot_fil_split, function(x) {
  suppressWarnings(
    x %>% 
      dplyr::select(-plot_cluster) %>%
      dplyr::select(where(~ any(. != 0))) %>%
      vegdist() %>%
      mean())
  }))

pdf(file = "img/plot_dist_hist.pdf", width = 8, height = 6)
ggplot() + 
  geom_histogram(data = data.frame(plot_dist_mean), aes(x = plot_dist_mean), 
    fill = pal[6], colour = "black", binwidth = 0.05) +
  geom_vline(xintercept = plot_dist_all_mean, colour = "red") + 
  theme_panel() + 
  labs(x = "Mean pairwise Bray distance within cluster", y = "Frequency")
dev.off()

# Exclude sites where Bray distance of plots is greater than mean of all pairs 
plot_dist_over <- which(plot_dist_mean < plot_dist_all_mean)
plot_dist_per <- round(length(plot_dist_over) / length(plot_dist_mean) * 100, 1)
plots_div_fil <- plots_div[plots_div$plot_cluster %in% names(plot_dist_over),]

# How many sites
n_sites <- nrow(plots_div_fil)

# Find plot IDs from filtered data 
plot_id_fil <- plots %>% 
  filter(plot_cluster %in% plots_div_fil$plot_cluster) %>% 
  pull(plot_id_vec) %>% 
  unlist() 

# How many trees were identified to species
ab_plot_mat_fil <- ab_plot_mat[rownames(ab_plot_mat) %in% plot_id_fil,]
nsp <- sum(ab_plot_mat_fil[,!grepl("indet", names(ab_plot_mat_fil))])
perspid <- nsp / sum(ab_plot_mat_fil) * 100

# Export TeX variables
write(
  c(
    commandOutput(plot_dist_per, "plotDistPer"),
    commandOutput(n_sites, "nSites"),
    commandOutput(clust_optim, "nCluster"),
    commandOutput(round(sil_mean_best, 2), "silBest"),
    commandOutput(perspid, "perSpID")
    ),
  file = "out/diversity_vars.tex")

# Write file
saveRDS(plots_div_fil, "dat/div.rds") 
  
