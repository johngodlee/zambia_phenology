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
library(xtable)

source("functions.R")

# Import data
dat <- readRDS("dat/plots_try.rds") 
dat <- st_as_sf(dat)

plot_id_lookup <- readRDS("dat/plot_id_lookup.rds")

ba_clust_mat <- readRDS("dat/ba_clust_mat.rds")

# Filter basal area matrix to plots in `dat`
ba_clust_mat <- ba_clust_mat[row.names(ba_clust_mat) %in% dat$plot_cluster,]
ba_clust_mat <- ba_clust_mat[,colSums(ba_clust_mat) > 0]

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
  slice_max(order_by = ba, n = 10, with_ties = FALSE) %>%
  rename(species_ba = ba ) %>%
  mutate(id = row_number()) %>%
  left_join(., sp_tax[,c("user_supplied_name", "family", "subfamily")], 
    by = c("species" = "user_supplied_name")) %>%
  rename(species_family = family, species_subfamily = subfamily) %>%
  mutate(species_prop_cum = cumsum(species_prop)) %>%
  pivot_wider(names_from = id,
    values_from = c("species", "species_ba", "species_prop", 
      "species_family", "species_subfamily", "species_prop_cum")) 

stopifnot(nrow(ba_clust_mat) == length(unique(dom_sp$plot_cluster)))

ba_chosen_mat <- ba_clust_mat

# Construct distance matrix
ba_dist <- vegdist(ba_chosen_mat, method = "bray")

# Ward distance
beta_ward <- hclust(ba_dist, method = "ward.D2")

# For each number of clusters create dataframe of dominant species
nclust <- seq_len(25)
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

# Identify optimal number of clusters
clust_optim <- nclust[sil_mean == max(sil_mean, na.rm = TRUE) & !is.na(sil_mean)]
#clust_optim <- 9

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

# How many of top n dominant species make up percentages of basal area?
pdf(file = "img/basal_area_dom_hist.pdf", width = 10, height = 8)
dom_sp %>%
  dplyr::select(contains("_cum_")) %>%
  gather(key, value, -plot_cluster) %>%
  ggplot(., aes(x = value)) + 
  geom_histogram(colour = "black", fill = "grey") + 
  geom_vline(xintercept = 0.8, linetype = 2, colour = "red") + 
  facet_wrap(~key) + 
  theme_panel()
dev.off()

# Reverse some traits so all increase with conservativism
rev_cols <- c()

cons_acq <- dat %>%
  st_drop_geometry() %>%
  dplyr::select(
    plot_cluster, 
    starts_with(c("leaf_n_mass", "leaf_p_mass", "sla", "wd"))) %>%
  mutate(across(all_of(rev_cols), 
      ~-1 * .x,
      .names = "{.col}_rev")) %>%
  dplyr::select(-all_of(rev_cols)) %>%
  filter(across(everything(), ~!is.na(.x)))

cons_pca <- prcomp(cons_acq[,-which(names(cons_acq) == "plot_cluster")], 
  center = TRUE, scale. = TRUE)

# Extract PCA values from each PCA
cons_pca_tidy <- as.data.frame(cons_pca$x) %>%
  mutate(plot_cluster = cons_acq$plot_cluster) %>%
  rename_with(~gsub("PC", "gen_pc", .x), starts_with("PC"))

cons_pca_arrows <- data.frame(x = rownames(cons_pca$rotation), 
  cons_pca$rotation)
cons_pca_arrows$x_lab <- case_when(
  cons_pca_arrows$x == "leaf_n_mass_genus_cwm" ~ "Leaf N",
  cons_pca_arrows$x == "leaf_p_mass_genus_cwm" ~ "Leaf P",
  cons_pca_arrows$x == "sla_genus_cwm" ~ "SLA",
  cons_pca_arrows$x == "wd_genus_cwm" ~ "WD",
  TRUE ~ NA_character_)

# Add values to data
dat_div <- dat %>%
  left_join(., div_df, by = "plot_cluster") %>%  # Diversity stats
  left_join(., dom_sp, by = "plot_cluster") %>%  # Dominant species
  left_join(., clust_best, by = "plot_cluster") %>%  # Ward clusters 
  left_join(., cons_pca_tidy, by = "plot_cluster") %>%  # Genus PCA cons-acq 
  rename(cluster = paste0("c", clust_optim))

pdf(file = "img/cons_pca_clust.pdf", width = 8, height = 5)
ggplot() + 
  geom_point(data = dat_div, 
    aes(x = gen_pc1, y = gen_pc2, fill = as.character(cluster)),
    shape = 21, colour = "black") + 
  stat_ellipse(data = dat_div,
    aes(x = gen_pc1, y = gen_pc2, colour = as.character(cluster)),
    show.legend = FALSE) + 
  geom_segment(data = cons_pca_arrows, aes(x = 0, y = 0, xend = (PC1*5),
    yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
    color = "black") +
  geom_label(data = cons_pca_arrows, aes(x = (PC1*4), y = (PC2*4),
    label = x_lab)) + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  theme_bw() + 
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  labs(x = "PC 1", y = "PC 2")
dev.off()

# Check no duplications or removals
stopifnot(nrow(dat) == nrow(dat_div))
stopifnot(all(!is.na(dat_div$plot_cluster)))

# Create table of species indicators and climatic data per cluster
clust_summ <- dat_div %>% 
  st_drop_geometry() %>%
  group_by(cluster) %>%
  summarise(
    richness_median = median(richness, na.rm = TRUE),
    richness_iqr = (quantile(richness, 0.75) - quantile(richness, 0.25)),
    n_sites = as.character(n()),
    map_mean = mean(map, na.rm = TRUE),
    map_sd = sd(map, na.rm = TRUE),
    diurnal_temp_range_mean = mean(diurnal_temp_range, na.rm = TRUE),
    diurnal_temp_range_sd = sd(diurnal_temp_range, na.rm = TRUE)) %>%
  mutate(
    cluster = as.character(cluster),
    richness = paste0(sprintf("%.0f", richness_median), "(", sprintf("%.0f", richness_iqr), ")"),
    map = paste0(sprintf("%.0f", map_mean), "(", sprintf("%.1f", map_sd), ")"),
    diurnal_temp_range = paste0(sprintf("%.0f", diurnal_temp_range_mean), "(", sprintf("%.1f", diurnal_temp_range_sd), ")")) %>%
  dplyr::select(cluster, n_sites, richness, map, diurnal_temp_range) %>%
  left_join(., indval_extrac_tidy, by = "cluster") %>%
  mutate(
    species = paste0("\\textit{", species, "}"),
    indval = sprintf("%.3f", indval))

clust_summ[ c(rbind(seq(1, nrow(clust_summ), 3), seq(3, nrow(clust_summ), 3))), 1:5] <- ""

names(clust_summ) <- c("Cluster", "N sites", "Richness", "MAP", "Diurnal $\\delta$T", "Species", "Indicator value")

# Export indval table
clust_summ_xtable <- xtable(clust_summ,
  label = "clust_summ",
  align = rep("c", 8),
  display = rep("s", 8),
  caption = "Climatic information and Dufrene-Legendre indicator species analysis for the vegetation type clusters identified by the PAM algorithm, based on basal area weighted species abundances. The three species per cluster with the highest indicator values are shown along with other key statistics for each cluster. MAP (Mean Annual Precipitation) and Diurnal $\\delta$T are reported as the mean and 1 standard deviation in parentheses. Species richness is reported as the median and the interquartile range in parentheses.")

fileConn <- file("out/clust_summ.tex")
writeLines(print(clust_summ_xtable, include.rownames = FALSE,
    table.placement = "H",
    hline.after = c(-1,0,seq(from = 3, by = 3, length.out = clust_optim)),
    sanitize.text.function = function(x) {x}),
  fileConn)
close(fileConn)

# Write file
saveRDS(dat_div, "dat/plots_div.rds") 
