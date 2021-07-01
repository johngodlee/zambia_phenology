# Get diversity statistics
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(tidyr)
library(sf)
library(vegan)
library(ade4)
library(ggplot2)
library(cluster)
library(labdsv)
library(shades)
library(tibble)
library(ggdendro)
library(taxize)

source("functions.R")

# Import data
plots <- readRDS("dat/plots_try.rds") 
plots <- st_as_sf(plots)

plot_id_lookup <- readRDS("dat/plot_id_lookup.rds")

ba_clust_mat <- readRDS("dat/ba_clust_mat.rds")

# Filter basal area matrix to plots in `plots`
ba_clust_mat <- ba_clust_mat[row.names(ba_clust_mat) %in% plots$plot_cluster,]
ba_clust_mat <- ba_clust_mat[,colSums(ba_clust_mat) > 0]

# Get taxonomic ranks
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

# Write species information file
saveRDS(sp_tax, "dat/sp_taxon.rds")

# Calculate common diversity statistics
div_df <- data.frame(plot_cluster = row.names(ba_clust_mat), 
  richness = unname(rowSums(ba_clust_mat != 0)),
  shannon = diversity(ba_clust_mat),
  simpson = diversity(ba_clust_mat, "simpson"),
  evenness = diversity(ba_clust_mat) / log(rowSums(ba_clust_mat > 0)))

# Effective true numbers diversity
# based on Shannon calculated by weighting on basal area.
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

# For each plot_cluster, dominant species by basal area
dom_sp <- do.call(rbind, lapply(ba_split, function(x) {
  ba_total <- sum(x$ba, na.rm = TRUE)
  x_order <- x[order(x$ba, decreasing = TRUE),]
  x_fil <- x_order[1:5,]
  x_fil$prop <- x_fil$ba / ba_total
  x_fil$id <- 1:5
  x_spread <- pivot_wider(x_fil, names_from = id, 
    values_from = c("species", "ba", "prop")) %>%
    filter(!is.na(plot_cluster))
  x_spread$ba_total <- ba_total

  return(x_spread)
}))

# Join families
dom_sp$family_1 <- sp_tax[match(dom_sp$species_1, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$family_2 <- sp_tax[match(dom_sp$species_2, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$family_3 <- sp_tax[match(dom_sp$species_3, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$family_4 <- sp_tax[match(dom_sp$species_4, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$family_5 <- sp_tax[match(dom_sp$species_5, sp_tax$user_supplied_name), "family"][[1]]

# Join subfamilies
dom_sp$subfamily_1 <- sp_tax[match(dom_sp$species_1, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$subfamily_2 <- sp_tax[match(dom_sp$species_2, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$subfamily_3 <- sp_tax[match(dom_sp$species_3, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$subfamily_4 <- sp_tax[match(dom_sp$species_4, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$subfamily_5 <- sp_tax[match(dom_sp$species_5, sp_tax$user_supplied_name), "subfamily"][[1]]

# Calculate cumulative proportions
dom_sp$prop_1_cum <- dom_sp$prop_1
dom_sp$prop_2_cum <- dom_sp$prop_1 + dom_sp$prop_2
dom_sp$prop_3_cum <- dom_sp$prop_2 + dom_sp$prop_3
dom_sp$prop_4_cum <- dom_sp$prop_3 + dom_sp$prop_4
dom_sp$prop_5_cum <- dom_sp$prop_4 + dom_sp$prop_5

# Gather dominant species dataframe
dom_gather <- dom_sp %>% 
  dplyr::select(plot_cluster, starts_with("species_")) %>%
  gather(dom, species, -plot_cluster, -ends_with("_cum")) %>%
  arrange(plot_cluster) %>%
  dplyr::select(-dom)

# Subset abundance matrix to only include dominant species per plot
ba_dom_mat <- ba_gather %>%
  inner_join(., dom_gather, by = c("plot_cluster", "species")) %>%
  spread(species, ba) %>%
  column_to_rownames("plot_cluster") %>%
  mutate(across(everything(), ~if_else(is.na(.x), 0, .x)))

# Check that dataframes are same size
stopifnot( nrow(ba_dom_mat[rowSums(ba_dom_mat) != 0,]) == 
  nrow(ba_clust_mat) )
stopifnot( ncol(dim(ba_dom_mat[,colSums(ba_dom_mat, na.rm = TRUE) != 0])) == 
  ncol(ba_clust_mat) )

# Construct distance matrix
ba_dist <- vegan::vegdist(ba_clust_mat)

# Ward distance
beta_ward <- hclust(ba_dist, method = "ward.D2")

# For each number of clusters create dataframe of dominant species
nclust <- seq_len(25)
dom_all <- lapply(nclust, function(y) {

  cuts <- cutree(beta_ward, y)

  sil <- silhouette(cuts, ba_dist)

  grp <- data.frame(
    plot_cluster = row.names(ba_clust_mat),
    cu = cuts)
  names(grp)[2] <- paste0("c", y)

  tmp <- aggregate(ba_clust_mat, list(grp[,2]), FUN = mean)
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

sil_mean <- unlist(lapply(sil_out, function(x) { 
  if (is.matrix(x)) {
    mean(x[,3], na.rm = TRUE)
  } else { 
    NA_real_
  }
}))

# Create silhouette plot
pdf(file = "img/ward_sil.pdf", width = 6, height = 4)
ggplot() + 
  geom_point(aes(y = sil_mean, x = nclust)) + 
  geom_line(aes(y = sil_mean, x = nclust)) + 
  theme_bw() + 
  labs(x = "N clusters", y = "Mean silhouette width")
dev.off()

# Create dataframe of dominant species
dom_out <- do.call(rbind, lapply(dom_all, "[[", 2))
clust_out <- Reduce( function(...) { merge(..., by = "plot_cluster") }, 
  lapply(dom_all, "[[", 1) )

# NSCA on basal area data
naxes <- 2
nsca <- dudi.nsc(df = ba_clust_mat, scannf = FALSE, nf = naxes)

nsca_inertia <- round(nsca$eig[naxes+1], 2) * 10

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
pdf(file = "img/nsca_sil.pdf", width = 6, height = 4)
ggplot(sil, aes(x = kval, y = v)) + 
  geom_point() + 
  geom_path() + 
  geom_vline(xintercept = n_clusters, linetype = 2) +
  theme_panel() + 
  labs(x = "N clusters", y = "Mean silhouette width") 
dev.off()

# Cluster euclidean distances of NSCA with ward algorithm
ward_clust <- hclust(nsca_dist, method = "ward.D2")
ward_dat <- dendro_data(ward_clust)

plot.new()
clust_classif <- rect.hclust(ward_clust, k = n_clusters)
nsca_df <- do.call(rbind, lapply(seq(length(clust_classif)), function(x) {
  data.frame(plot_cluster = names(clust_classif[[x]]), cluster = as.character(x))
}))

ward_merge <- left_join(ward_dat$labels, nsca_df, by = c("label" = "plot_cluster"))

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

# Indicator species per cluster
clust_indval <- indval(ba_clust_mat, clustering = nsc_df$cluster)

# Summarise indicator analysis
summary(clust_indval, p = 0.05, type = "short", digits = 2, show = p)
clust_indval$indval$sp <- row.names(clust_indval$indval)
row.names(clust_indval$indval) <- seq(from = 1, to = length(clust_indval$indval$sp))

nsc_indval_extrac <- lapply(1:n_clusters, function(x) {
    out <- head(clust_indval$indval[order(clust_indval$indval[[x]], 
          decreasing = TRUE),c(n_clusters+1, x)])
    out[!grepl("indet", out$sp, ignore.case = TRUE),]
  })

nsc_indval_extrac_tidy <- do.call(rbind, lapply(nsc_indval_extrac, function(x) {
    cluster <- names(x)[2]
    out <- x[1:3,]
    out$cluster <- cluster
    names(out) <- c("species", "indval", "cluster")
    out[,c(3,1,2)]
  })
)

# Combine clusters 1 and 3 based on indicator values?

# NSCA ordination plot with clusters
nsc_hull <- findHull(nsc_df, "RS1", "RS2", group = "cluster")

pdf(file = "img/nsca.pdf", width = 12, height = 6)
ggplot() +
  geom_polygon(data = nsc_hull, 
    aes(x = RS1, y = RS2, colour = cluster), fill = NA) + 
  geom_polygon(data = nsc_hull, 
    aes(x = RS1, y = RS2, fill = cluster), colour = NA, alpha = 0.5) + 
  geom_point(data = nsc_df, aes(x = RS1, y = RS2, fill = cluster), 
    shape = 21, size = 2) +
  scale_fill_manual(name = "Cluster", values = clust_pal) +
  scale_colour_manual(name = "Cluster", values = brightness(clust_pal, 0.5),
    limits = unique(nsc_df$cluster)) + 
  theme_panel() +
  labs(x = "NSC 1", y = "NSC 2")
dev.off()


# How many of top n dominant species make up percentages of basal area?
pdf(file = "img/basal_area_dom_hist.pdf", width = 10, height = 8)
dom_sp %>%
  dplyr::select(ends_with("_cum")) %>%
  gather(key, value, -plot_cluster) %>%
  ggplot(., aes(x = value)) + 
  geom_histogram(colour = "black", fill = "grey") + 
  geom_vline(xintercept = 0.8, linetype = 2, colour = "red") + 
  facet_wrap(~key) + 
  theme_panel()
dev.off()

# Reverse some traits so all increase with conservativism
rev_cols <- c("leaf_n_mass_cwm", "leaf_p_mass_cwm", "sla_cwm",
  "leaf_n_mass_genus_cwm", "leaf_p_mass_genus_cwm", "sla_genus_cwm")

cons_acq <- plots %>%
  st_drop_geometry() %>%
  dplyr::select(
    plot_cluster, 
    starts_with(c("leaf_n_mass", "leaf_p_mass", "sla", "ldmc", 
        "leaf_cn"))) %>%
  mutate(across(all_of(rev_cols), 
      ~-1 * .x,
      .names = "{.col}_rev")) %>%
  dplyr::select(-all_of(rev_cols)) %>%
  filter(across(everything(), ~!is.na(.x)))

# Get species and genus level conservative-acquisitive values
cons_acq_gen <- cons_acq %>%
  dplyr::select(plot_cluster, contains("_genus_cwm"))

cons_acq_sp <- cons_acq %>%
  dplyr::select(plot_cluster, !contains("_genus_cwm"))

# PCA of community weighted means of traits
cons_pca_sp <- prcomp(cons_acq_sp[,-which(names(cons_acq_sp) == "plot_cluster")], 
  center = TRUE, scale. = TRUE)

cons_pca_gen <- prcomp(cons_acq_gen[,-which(names(cons_acq_gen) == "plot_cluster")], 
  center = TRUE, scale. = TRUE)

# Extract PCA values from each PCA
cons_pca_sp_tidy <- as.data.frame(cons_pca_sp$x) %>%
  mutate(plot_cluster = cons_acq_sp$plot_cluster) %>%
  rename_with(~gsub("PC", "sp_pc", .x), starts_with("PC"))

cons_pca_sp_arrows <- data.frame(x = rownames(cons_pca_sp$rotation), 
  cons_pca_sp$rotation)

cons_pca_gen_tidy <- as.data.frame(cons_pca_gen$x) %>%
  mutate(plot_cluster = cons_acq_gen$plot_cluster) %>%
  rename_with(~gsub("PC", "gen_pc", .x), starts_with("PC"))

cons_pca_gen_arrows <- data.frame(x = rownames(cons_pca_gen$rotation), 
  cons_pca_gen$rotation)

cons_pca_sp_df <- dom_sp %>%
  left_join(., cons_pca_sp_tidy, by = "plot_cluster")

cons_pca_gen_df <- dom_sp %>%
  left_join(., cons_pca_gen_tidy, by = "plot_cluster") %>%
  mutate(genus_1 = unlist(lapply(strsplit(species_1, " "), "[", 1)))

# Plots of PCAs
# Species PCA
pdf(file = "img/cons_pca_sp.pdf", width = 20, height = 12)
ggplot() + 
  geom_point(data = cons_pca_sp_df, 
    aes(x = sp_pc1, y = sp_pc2, fill = species_1, size = prop_1, shape = species_1),
    colour = "black") +
  geom_segment(data = cons_pca_sp_arrows, aes(x = 0, y = 0, xend = (PC1*5),
    yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
    color = "black") +
  geom_label(data = cons_pca_sp_arrows, aes(x = (PC1*4), y = (PC2*4),
    label = x)) + 
  scale_shape_manual(values = rep(21:25, times = 14)) + 
  scale_fill_manual(values = rep(clust_pal, times = 14)) + 
  theme_bw() + 
  guides(fill = guide_legend(override.aes = list(size = 5))) 
dev.off()

# Genus PCA
pdf(file = "img/cons_pca_gen.pdf", width = 20, height = 12)
ggplot() + 
  geom_point(data = cons_pca_gen_df, 
    aes(x = gen_pc1, y = gen_pc2, fill = genus_1, size = prop_1, shape = genus_1), 
    colour = "black") + 
  geom_segment(data = cons_pca_gen_arrows, aes(x = 0, y = 0, xend = (PC1*5),
    yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
    color = "black") +
  geom_label(data = cons_pca_gen_arrows, aes(x = (PC1*4), y = (PC2*4),
    label = x)) + 
  scale_shape_manual(values = rep(21:25, times = 14)) + 
  scale_fill_manual(values = rep(clust_pal, times = 14)) + 
  theme_bw() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()

# Add values to data
plots_div <- plots %>%
  left_join(., div_df, by = "plot_cluster") %>%  # Diversity stats
  left_join(., dom_sp, by = "plot_cluster") %>%  # Dominant species
  left_join(., nsc_df, by = "plot_cluster") %>%  # NSCA clusters and axes
  left_join(., cons_pca_sp_tidy, by = "plot_cluster") %>%  # Species PCA cons-acq
  left_join(., cons_pca_gen_tidy, by = "plot_cluster")  # Genus PCA cons-acq

# Check no duplications or removals
stopifnot(nrow(plots) == nrow(plots_div))
stopifnot(all(!is.na(plots_div$plot_cluster)))

# Write file
saveRDS(plots_div, "dat/plots_div.rds") 
