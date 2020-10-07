# Get diversity statistics
# John Godlee (johngodlee@gmail.com)
# 2020-09-02

# Packages
library(dplyr)
library(sf)
library(vegan)
library(ape)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(viridis)
library(ade4)
library(cluster)
library(labdsv)
library(shades)

source("functions.R")

# Import data
tree_ab_mat <- readRDS("dat/tree_ab_mat.rds")
tree_ab_mat_plot <- readRDS("dat/tree_ab_mat_plot.rds")

plots <- readRDS("dat/plots_phen.rds") 

plot_id_lookup <- readRDS("dat/plot_id_lookup.rds")

# Exclude clusters with mental species composition
# Exclude species with only one occurrence
ab_mat_clean <- tree_ab_mat %>%
  filter(!row.names(.) %in% c("ZIS_2385", "ZIS_3765"),
    row.names(.) %in% unique(plots$plot_cluster)) %>%
  dplyr::select(which(!colSums(., na.rm = TRUE) %in% c(0,1)))

# Conduct NSCA (Non-symmetric Correspondence Analysis) 
##' 4 axes
nsca <- dudi.nsc(df = ab_mat_clean, scannf = FALSE, nf = 4)
##' p = species
##' n = sites

# Run clustering around medoids
##' 4 clusters
pam_clust <- pam(nsca$l1, k = 4, metric = "manhattan", diss = FALSE)

# Extract data from clustering and NSCA
nsc_df <- data.frame(plot_cluster = names(pam_clust$clustering), 
  cluster = unname(pam_clust$clustering)) %>% 
  left_join(., rownames_to_column(nsca$li), by = c("plot_cluster" = "rowname")) %>%
  mutate(cluster = as.character(cluster)) %>%
  rename_with(~gsub("Axis", "_", paste0("nsc", .)), starts_with("Axis"))

# Plot of clusters
nscaPlot <- function(x,y, dat = nsc_df) {
  x <- ensym(x)
  y <- ensym(y)

  p <- ggplot(dat, aes(x = !!x, y = !!y)) + 
    geom_point(shape = 21, size = 2, aes(fill = cluster)) + 
    geom_bag(prop = 0.95, alpha = 0.2, aes(colour = cluster, fill = cluster)) + 
    scale_fill_manual(name = "Cluster", values = clust_pal) + 
    scale_x_continuous(breaks = seq(-10,10, by = 2)) +
    scale_y_continuous(breaks = seq(-10,10, by = 2)) +
    theme_panel() + 
    labs(x = gsub("nsc_", "NSC ", x),
      y = gsub("nsc_", "NSC ", y))
}

pdf(file = "img/nsca.pdf", width = 12, height = 6)
nscaPlot(nsc_1, nsc_2) + 
  nscaPlot(nsc_3, nsc_4) + 
  plot_layout(guides = "collect", widths = 1) & 
  scale_colour_manual(name = "Cluster", values = brightness(clust_pal, 0.5), 
    limits = unique(nsc_df$cluster))
dev.off()

# Indicator species per cluster
clust_indval <- indval(ab_mat_clean, clustering = nsc_df$cluster)

# Summarise indicator analysis
summary(clust_indval, p = 0.05, type = "short", digits = 2, show = p)
clust_indval$indval$sp <- row.names(clust_indval$indval)
row.names(clust_indval$indval) <- seq(from = 1, to = length(clust_indval$indval$sp))

indval_extrac <- lapply(1:4, function(x) {
    out <- head(clust_indval$indval[order(clust_indval$indval[[x]], 
          decreasing = TRUE),c(5, x)])
    out[!grepl("indet", out$sp, ignore.case = TRUE),]
  })

indval_extrac_tidy <- do.call(rbind, lapply(indval_extrac, function(x) {
    cluster <- names(x)[2]
    out <- x[1:3,]
    out$cluster <- cluster
    names(out) <- c("Species", "Indicator value", "Cluster")
    out[,c(3,1,2)]
  })
)

# Export indval table
indval_xtable <- xtable(indval_extrac_tidy, 
  label = "indval",
  digits = 3,
  align = c("c", "c", "c", "c"),
  display = c("s", "s", "s", "f"))

fileConn <- file("out/indval.tex")
writeLines(print(indval_xtable, include.rownames = FALSE,
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Calculate common  diversity statistics
div_df <- data.frame(plot_cluster = row.names(ab_mat_clean), 
  richness = unname(rowSums(ab_mat_clean != 0)),
  shannon = diversity(ab_mat_clean),
  simpson = diversity(ab_mat_clean, "simpson"),
  evenness = diversity(ab_mat_clean, "invsimpson") / rowSums(ab_mat_clean > 0))

# Add values to data
plots_div <- plots %>%
  left_join(., nsc_df, by = "plot_cluster") %>%
  left_join(., div_df, by = "plot_cluster") %>% 
  filter(!is.na(richness))

saveRDS(plots_div, "dat/plots_div.rds") 

# Test Beta diversity between plots in a cluster - UNFINISHED

# Filter to plots we're using
tree_ab_mat_plot_fil <- tree_ab_mat_plot[row.names(tree_ab_mat_plot) %in% 
  unlist(strsplit(plots_div$plot_id, ",")) & rowSums(tree_ab_mat_plot) != 0,]

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
plot_dist_per <- length(which(plot_dist_mean_clean < plot_dist_all_mean)) / 
  length(plot_dist_mean_clean) * 100

write(commandOutput(plot_dist_per, "plotDistPer"), file="out/vars.tex", 
  append = TRUE)

