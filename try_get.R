# Extract TRY species and trait codes
# John Godlee (johngodlee@gmail.com)
# 2021-02-05

# Packages
library(dplyr)
library(tidyr)
library(data.table)
library(BIOMASS)
library(fst)

# Import data
plots <- readRDS("dat/plots.rds")
ba_mat <- readRDS("dat/ba_clust_mat.rds")
try_traits <- read.table("dat/try/tde202125164358.txt", 
  skip = 3, header = TRUE, sep = "\t")
try_sp <- read.table("dat/try/TryAccSpecies.txt",
  header = TRUE, sep = "\t")
# try_dat <- fread("/Volumes/seosaw_spat/try/14425.txt", 
#   header = TRUE, sep = "\t", dec = ".", quote = "", data.table = FALSE)
# write.fst("/Volumes/seosaw_spat/try/14425.fst")
try_dat <- read.fst("/Volumes/seosaw_spat/try/14425.fst")

# Extract genera from abundance matrix
genus <- unique(unlist(lapply(strsplit(names(ba_mat), " "), "[", 1)))

# Create genus column in TRY species list
try_sp$genus <- unlist(lapply(strsplit(try_sp$AccSpeciesName, " "), "[", 1))

# Find all species from genera in the plot data
genus_out <- try_sp[try_sp$genus %in% genus,
  c("AccSpeciesName", "AccSpeciesID")]
names(genus_out) <- c("species", "id")

# Remove genus only names
all_sp_out <- genus_out[
  unlist(lapply(strsplit(genus_out$species, " "), length)) > 1,]

# Find traits
trait_out <- try_traits %>%
  filter(TraitID %in% c(47,3115,3116,3117,14,46,8,50,15,146,1229,24)) %>%
  arrange(desc(ObsNum)) %>%
  pull(TraitID) %>%
  as.character()

# Write files
write.csv(all_sp_out, "dat/try/species_list.csv", row.names = FALSE)
writeLines(trait_out, "dat/try/trait_list.txt")

# Remove genera only from species list

# Subset TRY data to species and columns of interest
try_clean <- try_dat %>%
  dplyr::select(
    dataset = Dataset,
    obs_id = ObservationID,
    species_id = AccSpeciesID,
    species = AccSpeciesName,
    trait_id = TraitID,
    trait = TraitName,
    key_id = DataID,
    key = DataName,
    val = OrigValueStr,
    unit = OrigUnitStr,
    val_std = StdValue,
    unit_std = UnitName,
    error_risk = ErrorRisk) %>%
  filter(species_id %in% all_sp_out$id)

# Trait ID match
trait_id_lookup <- try_clean %>%
  dplyr::select(trait_id, trait) %>%
  filter(
    !is.na(trait_id),
    trait_id %in% trait_out) %>%
  distinct() %>%
  mutate(
    trait = case_when(
      trait_id %in% c(3115, 3116, 3117) ~ gsub(":.*", "", trait),
      TRUE ~ trait),
    trait_short = case_when(
      trait_id == 8 ~ "N fixation capacity",
      trait_id == 14 ~ "Leaf N / dry mass",
      trait_id == 15 ~ "Leaf P / dry mass",
      trait_id == 24 ~ "Bark thickness",
      trait_id == 46 ~ "Leaf thickness",
      trait_id == 47 ~ "Leaf dry matter content",
      trait_id == 50 ~ "Leaf N / leaf area",
      trait_id == 146 ~ "Leaf C:N",
      trait_id == 1229 ~ "Wood N / dry mass",
      trait_id %in% c(3115, 3116, 3117) ~ "SLA",
      TRUE ~ as.character(trait_id)),
    trait_col = case_when(
      trait_id == 8 ~ "n_fix_cap",
      trait_id == 14 ~ "leaf_n_mass",
      trait_id == 15 ~ "leaf_p_mass",
      trait_id == 24 ~ "bark_thick",
      trait_id == 46 ~ "leaf_thick",
      trait_id == 47 ~ "ldmc",
      trait_id == 50 ~ "leaf_n_area",
      trait_id == 146 ~ "leaf_cn",
      trait_id == 1229 ~ "wood_n_mass",
      trait_id %in% c(3115, 3116, 3117) ~ "sla",
      TRUE ~ as.character(trait_id))
  )

write.csv("dat/trait_id_lookup.csv", row.names = FALSE)

# Split try data by observation
try_split <- split(try_clean, try_clean$obs_id)

# Create clean dataframe of observations
try_species <- as.data.frame(do.call(rbind, lapply(try_split, function(x) {
  if(any(x$trait_id %in% trait_id_lookup$trait_id)) {
    traits <- x[x$trait_id %in% trait_id_lookup$trait_id,
      c("species", "trait_id", "val_std")]

    traits$trait_id <- trait_id_lookup$trait_col[
    match(traits$trait_id, trait_id_lookup$trait_id)]

    return(traits)
  } else {
    return(NULL)
  }
})))
names(try_species) <- c("binom", "trait_id", "val_std")

# Get wood density
wd_sp <- wdData %>%
  mutate(binom = paste(genus, species)) %>%
  filter(
    genus %in% genus,
    grepl("africa", region, ignore.case = TRUE)) %>%
  group_by(binom) %>%
  summarise(val_std = mean(wd, na.rm = TRUE)) %>%
  mutate(trait_id = "wd")

# Bind all traits for all species in the same genera as in plot data
traits_sp <- bind_rows(wd_sp, try_species) %>%
  filter(!is.na(val_std)) %>%
  group_by(binom, trait_id) %>%
  summarise(val_std = mean(val_std, na.rm = TRUE)) %>%
  spread(trait_id, val_std)


# Write to file
saveRDS(traits_sp, "dat/traits_sp.rds")
