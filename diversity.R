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

# Write species information to file
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

# Gather basal area matrix to dataframe
ba_gather <- ba_clust_mat %>%
  rownames_to_column("plot_cluster") %>%
  gather(species, ba, -plot_cluster) %>%
  filter(ba > 0) %>%
  group_by(plot_cluster, species) %>%
  summarise(ba = sum(ba, na.rm = TRUE)) 

# Split by plot_cluster
ba_split <- split(ba_gather, ba_gather$plot_cluster)

# For each plot_cluster, dominant species and genera by basal area
dom_sp <- ba_gather %>% 
  group_by(plot_cluster) %>%
  mutate(species_prop = ba / sum(ba, na.rm = TRUE)) %>%
  slice_max(order_by = ba, n = 5) %>%
  rename(species_ba = ba ) %>%
  mutate(id = row_number()) %>%
  pivot_wider(names_from = id,
    values_from = c("species", "species_ba", "species_prop")) 

dom_sp$species_family_1 <- sp_tax[match(dom_sp$species_1, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$species_family_2 <- sp_tax[match(dom_sp$species_2, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$species_family_3 <- sp_tax[match(dom_sp$species_3, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$species_family_4 <- sp_tax[match(dom_sp$species_4, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$species_family_5 <- sp_tax[match(dom_sp$species_5, sp_tax$user_supplied_name), "family"][[1]]
dom_sp$species_subfamily_1 <- sp_tax[match(dom_sp$species_1, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$species_subfamily_2 <- sp_tax[match(dom_sp$species_2, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$species_subfamily_3 <- sp_tax[match(dom_sp$species_3, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$species_subfamily_4 <- sp_tax[match(dom_sp$species_4, sp_tax$user_supplied_name), "subfamily"][[1]]
dom_sp$species_subfamily_5 <- sp_tax[match(dom_sp$species_5, sp_tax$user_supplied_name), "subfamily"][[1]]

dom_sp$species_prop_1_cum <- dom_sp$species_prop_1
dom_sp$species_prop_2_cum <- dom_sp$species_prop_1 + dom_sp$species_prop_2
dom_sp$species_prop_3_cum <- dom_sp$species_prop_2 + dom_sp$species_prop_3
dom_sp$species_prop_4_cum <- dom_sp$species_prop_3 + dom_sp$species_prop_4
dom_sp$species_prop_5_cum <- dom_sp$species_prop_4 + dom_sp$species_prop_5

dom_gen <- ba_gather %>% 
  mutate(genus = gsub(" .*", "", species)) %>%
  group_by(plot_cluster, genus) %>%
  summarise(genus_ba = sum(ba, na.rm = TRUE)) %>%
  group_by(plot_cluster) %>%
  mutate(genus_prop = genus_ba / sum(genus_ba, na.rm = TRUE)) %>%
  slice_max(order_by = genus_ba, n = 5) %>%
  mutate(id = row_number()) %>%
  pivot_wider(names_from = id,
    values_from = c("genus", "genus_ba", "genus_prop")) 

dom_sp_gen <- full_join(dom_sp, dom_gen, "plot_cluster")

# Abundance matrix by dominant genera only
dom_gen_gather <- dom_gen %>%
  dplyr::select(1:6) %>%
  gather(key, genus, -plot_cluster) %>% 
  dplyr::select(-key)

ba_dom_mat <- ba_gather %>%
  mutate(genus = gsub(" .*", "", species)) %>%
  dplyr::select(-species) %>%
  group_by(plot_cluster, genus) %>%
  mutate(ba = sum(ba, na.rm = TRUE)) %>%
  inner_join(., dom_gen_gather, by = c("plot_cluster", "genus")) %>%
  distinct() %>%
  spread(genus, ba) %>%
  column_to_rownames("plot_cluster") %>%
  mutate(across(everything(), ~if_else(is.na(.x), 0, .x))) %>% 
  dplyr::select(which(colSums(.) > 0))

# Construct distance matrix
ba_dist <- vegan::vegdist(ba_dom_mat)

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

clust_optim <- nclust[sil_mean == max(sil_mean, na.rm = TRUE) & !is.na(sil_mean)]
clust_optim <- 4

# Create silhouette plot
pdf(file = "img/ward_sil.pdf", width = 6, height = 4)
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
clust_indval <- indval(ba_dom_mat, clustering = clust_best[[paste0("c", clust_optim)]])

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

cons_pca_gen <- prcomp(cons_acq_gen[,-which(names(cons_acq_gen) == "plot_cluster")], 
  center = TRUE, scale. = TRUE)

# Extract PCA values from each PCA
cons_pca_gen_tidy <- as.data.frame(cons_pca_gen$x) %>%
  mutate(plot_cluster = cons_acq_gen$plot_cluster) %>%
  rename_with(~gsub("PC", "gen_pc", .x), starts_with("PC"))

cons_pca_gen_arrows <- data.frame(x = rownames(cons_pca_gen$rotation), 
  cons_pca_gen$rotation)

cons_pca_gen_df <- dom_gen %>%
  left_join(., cons_pca_gen_tidy, by = "plot_cluster")

# Genus PCA plot
pdf(file = "img/cons_pca_gen.pdf", width = 20, height = 12)
ggplot() + 
  geom_point(data = cons_pca_gen_df, 
    aes(x = gen_pc1, y = gen_pc2, fill = genus_1, size = genus_prop_1, shape = genus_1), 
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
  left_join(., clust_best, by = "plot_cluster") %>%  # Ward clusters 
  left_join(., cons_pca_gen_tidy, by = "plot_cluster") %>%  # Genus PCA cons-acq 
  rename(cluster = paste0("c", clust_optim))

pdf(file = "img/cons_pca_clust.pdf", width = 8, height = 5)
ggplot() + 
  geom_point(data = plots_div, 
    aes(x = gen_pc1, y = gen_pc2, 
      colour = as.character(cluster))) + 
  geom_segment(data = cons_pca_gen_arrows, aes(x = 0, y = 0, xend = (PC1*5),
    yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
    color = "black") +
  geom_label(data = cons_pca_gen_arrows, aes(x = (PC1*4), y = (PC2*4),
    label = x)) + 
  theme_bw() + 
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()


# Check no duplications or removals
stopifnot(nrow(plots) == nrow(plots_div))
stopifnot(all(!is.na(plots_div$plot_cluster)))

# Write file
saveRDS(plots_div, "dat/plots_div.rds") 
