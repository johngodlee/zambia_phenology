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
library(viridis)

source("functions.R")

# Import data
tree_ab_mat <- readRDS("dat/tree_ab_mat.rds")

plots <- readRDS("dat/plots_phen.rds") 

# Exclude plot with mental species composition
# Exclude species with only one occurrence
ab_mat_clean <- tree_ab_mat %>%
  dplyr::select(which(!colSums(., na.rm=TRUE) %in% 1)) %>%
  filter(!row.names(.) %in% c("ZIS_2385", "ZIS_3765"),
    row.names(.) %in% unique(plots$plot_cluster)) 

# NMDS dimensions
pdf(file = "img/nmds_scree.pdf", width = 5, height = 5)
NMDS.scree(ab_mat_clean, trys = 10, distance = "jaccard")
dev.off()

# Run NMDS
##' 9 dimensions keeps stress below 0.1
##' 3 dimensions keeps below 0.15
nmds <- metaMDS(ab_mat_clean,
  autotransform = TRUE, distance = "bray", try = 20, trymax = 20, 
  k = 4, stepacross = TRUE)

# Check output
pdf(file = "img/nmds_stress.pdf", width = 5, height = 5)
stressplot(nmds)
dev.off()

# Extract Plot scores from NMDS
plot_scores <- as.data.frame(scores(nmds))  
plot_scores$plot_cluster <- row.names(plot_scores)
row.names(plot_scores) <- NULL

# Extract species scores from NMDS analysis
species_scores <- as.data.frame(scores(nmds, "species")) 
species_scores$species_binomial <- rownames(species_scores)
row.names(species_scores) <- NULL

# Save NMDS 
saveRDS(species_scores, "dat/species_scores.rds")
saveRDS(nmds, "dat/nmds.rds")

# Calculate common  diversity statistics
div_df <- data.frame(plot_cluster = row.names(tree_ab_mat), 
  richness = unname(rowSums(tree_ab_mat != 0)),
  shannon = diversity(tree_ab_mat),
  simpson = diversity(tree_ab_mat, "simpson"),
  evenness = diversity(tree_ab_mat, "invsimpson") / rowSums(tree_ab_mat > 0))

# Add values to data
plots_div <- plots %>%
  left_join(., plot_scores, by = "plot_cluster") %>%
  left_join(., div_df, by = "plot_cluster") %>% 
  filter(!is.na(richness))

saveRDS(plots_div, "dat/plots_div.rds") 

