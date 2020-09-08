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

ab_mat_clean <- tree_ab_mat[row.names(tree_ab_mat) != "ZIS_2385",]

plots <- readRDS("dat/plots.rds") 

# NMDS
NMDS.scree <- function(x, dims = 10, ...) {
  # Create dimensions vector
  if (length(dims) == 1) {
    dim_vec <- seq(dims) 
  } else {
    dim_vec <- dims
  }

  # Create list of metaMDS objects
  meta_list <- lapply(dim_vec, function(y) {
    metaMDS(x, autotransform = F, k = y, ...)
  })

  # Extract stress values
  stress_vec <- unname(unlist(lapply(meta_list, `[`, "stress")))

  # Create plot
  plot(dim_vec, stress_vec,
    xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
}

NMDS.scree(ab_mat_clean, trys = 10)

tree_ab_nmds <- metaMDS(ab_mat_clean, distance = "bray", try = 100, 
  trymax = 100, k = 4, autotransform = FALSE)

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

# Traits TESTING!!!!
library(BIEN)
library(purrr)

# Select traits
traits_all <- BIEN_trait_list()
traits_sel <- traits_all[c(6,8,16,24,43,49),]

# Species names
species <- names(tree_ab_mat)
genera_split <- strsplit(names(tree_ab_mat), " ")
genera <- unique(sapply(genera_split, `[`, 1))

trait_species <- BIEN_trait_traitbyspecies(species = species, trait = traits_sel)
trait_genus <- BIEN_trait_traitbygenus(genus = genera, trait = traits_sel)

trait_species_split <- split(trait_species, trait_species$trait_name)
trait_genus_split <- split(trait_genus, trait_genus$trait_name)

bien_clean <- function(x, name) {
  traits_summ <- x %>%
    group_by(scrubbed_species_binomial) %>%
    summarise(
      trait_name = first(na.omit(trait_name)),
      trait_value = mean(as.numeric(trait_value), na.rm = TRUE),
      unit = first(na.omit(unit))) %>%
    mutate(trait_name = paste(trait_name, unit, sep = "_")) %>%
    dplyr::select(-unit)
  colname <- gsub("\\s", "_", 
    gsub("\\.", "_", 
      gsub("\\-", "", unique(traits_summ$trait_name))))
  traits_clean <- traits_summ[,c(1,3)]
  names(traits_clean) <- c(name, colname)
  return(traits_clean)
  }

traits_species_summ <- lapply(trait_species_split[1:6], bien_clean, name = "species")
traits_genus_summ <- lapply(trait_genus_split[1:6], bien_clean, name = "genus")
traits_df <- reduce(traits_summ, left_join, by = "species") %>%
  right_join(., data.frame(species), by = "species")

traits_df$genus <- unique(sapply(strsplit(traits_df$species, " "), `[`, 1))
traits_df$species <- unique(sapply(strsplit(traits_df$species, " "), `[`, 2))

traits_genus_summ


# How many aren't NA?
apply(traits_df, 2, function(x) {
  length(which(!is.na(x)))
  })
