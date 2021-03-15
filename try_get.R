# Extract TRY species and trait codes
# John Godlee (johngodlee@gmail.com)
# 2021-02-05

# Packages
library(dplyr)
library(data.table)

# Import data
ab <- readRDS("dat/ba_clust_mat.rds")
try_traits <- read.table("dat/try/tde202125164358.txt", 
  skip = 3, header = TRUE, sep = "\t")
try_sp <- read.table("dat/try/TryAccSpecies.txt",
  header = TRUE, sep = "\t")

# Extract genera from abundance matrix
genus <- unique(unlist(lapply(strsplit(names(ab), " "), "[", 1)))

# Create genus column in TRY species list
try_sp$genus <- unlist(lapply(strsplit(try_sp$AccSpeciesName, " "), "[", 1))

# Match genus names
genus_out <- try_sp[try_sp$genus %in% genus,
  c("AccSpeciesName", "AccSpeciesID")]
names(genus_out) <- c("species", "id")

# Find traits
trait_out <- try_traits %>%
  filter(TraitID %in% c(47,3115,3116,3117,14,46,8,50,15,146,1229,24)) %>%
  arrange(desc(ObsNum)) %>%
  pull(TraitID) %>%
  as.character()

# Write files
write.csv(genus_out, "dat/try/species_list.csv", row.names = FALSE)
writeLines(trait_out, "dat/try/trait_list.txt")
