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

plots <- readRDS("dat/plots.rds") 

# Exclude plot with mental species composition
t(tree_ab_mat[row.names(tree_ab_mat) == "ZIS_2385",]) %>%
  as.data.frame() %>%
  filter(ZIS_2385 > 0)
##" Only Pterocarpus lucens
ab_mat_clean <- tree_ab_mat[row.names(tree_ab_mat) != "ZIS_2385",]

# NMDS dimensions
pdf(file = "img/nmds_scree.pdf", width = 5, height = 5)
NMDS.scree(ab_mat_clean, trys = 10)
dev.off()

# Run NMDS
nmds <- metaMDS(ab_mat_clean, distance = "jaccard", try = 100, 
  trymax = 150, k = 4, autotransform = FALSE)

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

