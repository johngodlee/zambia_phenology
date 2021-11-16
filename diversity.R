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

ab_plot_mat <- readRDS("dat/ab_plot_mat.rds")

ba_clust_mat <- readRDS("dat/ba_clust_mat.rds")

evi <- readRDS("dat/evi_pred_all.rds")

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

stopifnot(all(!duplicated(ba_perc_subfab_wide$plot_cluster)))

saveRDS(ba_perc_fam, "dat/ba_perc_fam.rds")
saveRDS(ba_perc_subfab_wide, "dat/ba_perc_subfab.rds")

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

sil_mean_best <- max(sil_mean, na.rm = TRUE) * 10

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
  left_join(., ba_perc_subfab_wide, by = "plot_cluster") %>%  # Dominant species
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

# Calculate mean pairwise Bray distances within each cluster
# Filter to plots we're using 
tree_ab_mat_plot_fil <- ab_plot_mat %>%
  rownames_to_column("plot_id") %>%
  left_join(., plot_id_lookup, by = "plot_id") %>%
  filter(plot_cluster %in% dat_div$plot_cluster) %>%
  dplyr::select(-plot_id) 

# Mean pairwise Bray distance across all plots
tree_ab_mat_plot_fil_fil <- tree_ab_mat_plot_fil %>%
  dplyr::select(-plot_cluster) %>%
  filter(rowSums(.) != 0)
plot_dist_all_mean <- mean(vegdist(tree_ab_mat_plot_fil_fil))

# Split by plot cluster
tree_ab_mat_plot_fil_split <- split(tree_ab_mat_plot_fil, tree_ab_mat_plot_fil$plot_cluster)

# Bray distance among plots in each cluster
plot_dist_mean <- unlist(lapply(tree_ab_mat_plot_fil_split, function(x) {
  x %>% 
    dplyr::select(-plot_cluster) %>%
    dplyr::select(where(~ any(. != 0))) %>%
    vegdist() %>%
    mean()
  }))

pdf(file = "img/plot_dist_hist.pdf", width = 8, height = 6)
ggplot() + 
  geom_histogram(data = data.frame(plot_dist_mean), aes(x = plot_dist_mean), 
    fill = pal[6], colour = "black", binwidth = 0.05 ) +
  geom_vline(xintercept = plot_dist_all_mean, colour = "red") + 
  theme_panel() + 
  labs(x = "Mean pairwise Bray distance within cluster", y = "Frequency")
dev.off()

# Percentage of clusters with mean pairwise distance lower than the mean across all pairs
plot_dist_mean_clean <- plot_dist_mean[!is.na(plot_dist_mean)]
plot_dist_over <- which(plot_dist_mean_clean < plot_dist_all_mean)
plot_dist_per <- round(length(plot_dist_over) / length(plot_dist_mean_clean) * 100, 1)

dat_div_fil <- dat_div[dat_div$plot_cluster %in% names(plot_dist_over),]

# How many sites
n_sites <- nrow(dat_div_fil)

# Create table of species indicators and climatic data per cluster
clust_summ <- dat_div_fil %>% 
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

names(clust_summ) <- c("Cluster", "N sites", "Richness", "MAP", "$\\delta$T", "Species", "Indicator value")

# Export indval table
clust_summ_xtable <- xtable(clust_summ,
  label = "clust_summ",
  align = rep("c", 8),
  display = rep("s", 8),
  caption = "Climatic information and Dufrene-Legendre indicator species analysis for the vegetation type clusters identified by the PAM algorithm, based on basal area weighted species abundances. The three species per cluster with the highest indicator values are shown along with other key statistics for each cluster. MAP (Mean Annual Precipitation) and $\\delta$T (Diurnal temperature range) are reported as the mean and 1 standard deviation in parentheses. Species richness is reported as the median and the interquartile range in parentheses.")

fileConn <- file("out/clust_summ.tex")
writeLines(print(clust_summ_xtable, include.rownames = FALSE,
    table.placement = "H",
    hline.after = c(-1,0,seq(from = 3, by = 3, length.out = clust_optim)),
    sanitize.text.function = function(x) {x}),
  fileConn)
close(fileConn)


# Subset EVI time series to plots 
evi_fil <- evi[names(evi) %in% dat_div_fil$plot_cluster]

stopifnot(all(dat_div_fil$plot_cluster == names(evi_fil)))

# Get all time series and clusters in tidy dataframe
gam_all <- do.call(rbind, lapply(unique(dat_div_fil$cluster), function(x) {
  out_df <- do.call(rbind, lapply(evi_fil[names(evi_fil) %in% 
      dat_div_fil$plot_cluster[dat_div_fil$cluster == x]], function(y) {
    out <- y$seas
    out$plot_cluster <- unique(y$series$plot_cluster)
    return(out)
  }))
  out_df$cluster <- x
  return(out_df)
}))

gam_all_clean <- gam_all[gam_all$doy < quantile(gam_all$doy, 0.99) &
  gam_all$doy > quantile(gam_all$doy, 0.01),]

# Plot of mean time series per cluster
pdf(file = "img/gam_compare_clust.pdf", width = 8, height = 5)
ggplot() + 
  stat_summary(data = gam_all_clean, 
     aes(x = doy, y = pred, fill = as.character(cluster)),
     fun.data = mean_se, geom = "ribbon", alpha = 0.5) + 
   stat_summary(data = gam_all_clean, 
     aes(x = doy, y = pred, colour = as.character(cluster)), 
     fun = mean, geom = "line") + 
  scale_fill_manual(name = "Cluster", values = clust_pal) +
  scale_colour_manual(name = "Cluster", values = clust_pal) +
  theme_panel() + 
  labs(x = "Days from 1st Jan.", y = "EVI") 
dev.off()

write(
  c(
    commandOutput(plot_dist_per, "plotDistPer"),
    commandOutput(n_sites, "nSites"),
    commandOutput(clust_optim, "nCluster"),
    commandOutput(round(sil_mean_best, 2), "silBest")
    ),
  file = "out/diversity_vars.tex")

# Write file
saveRDS(dat_div_fil, "dat/plots_div.rds") 
