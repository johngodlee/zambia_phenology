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

# Run PCOA on distance matrix
tree_dist <- vegdist(tree_ab_mat, method = "bray")
tree_pcoa <- pcoa(tree_dist)

# Percentage variance explained
eig_per <- list()
for (i in 1:length(tree_pcoa$values$Relative_eig)) {
  eig_per[i] <- sum(tree_pcoa$values$Relative_eig[1:i]) / sum(tree_pcoa$values$Relative_eig)
}

# First three axes
write(
  commandOutput(numFormat(sum(unlist(eig_per[1:3]))), "pcoaPer"),
  file="out/vars.tex", append=TRUE)

# Extract arrow vectors for species
trait_pcoa_arrows <- pcoaArrows(tree_pcoa, tree_ab_mat, sort = TRUE)

arrows_main <- trait_pcoa_arrows$U
arrows_main$species <- row.names(arrows_main)
arrows_main$species <- gsub(" ", "\n", arrows_main$species)

saveRDS(arrows_main, "dat/pcoa_arrows.rds")

# Extract values
pcoa_values <- as.data.frame(tree_pcoa$vectors)[,1:5]
names(pcoa_values) <- paste("pcoa", 1:5, sep = "_")
pcoa_values$plot_cluster <- row.names(pcoa_values)

# Calculate common  diversity statistics
div_df <- data.frame(plot_cluster = row.names(tree_ab_mat), 
  richness = unname(rowSums(tree_ab_mat != 0)),
  shannon = diversity(tree_ab_mat),
  simpson = diversity(tree_ab_mat, "simpson"),
  evenness = diversity(tree_ab_mat, "invsimpson") / rowSums(tree_ab_mat > 0))

# Add values to data
plots_div <- plots %>%
  left_join(., pcoa_values, by = "plot_cluster") %>%
  left_join(., div_df, by = "plot_cluster") %>% 
  filter(!is.na(richness))

saveRDS(plots_div, "dat/plots_div.rds") 

