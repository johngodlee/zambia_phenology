# Extract TRY database traits - conservative/acquisitive strategies
# John Godlee (johngodlee@gmail.com)
# 2021-02-06

# Packages
library(sf)
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)

# Import data
try_dat <- fread("dat/try/13583.txt", 
  header = TRUE, sep = "\t", dec = ".", quote = "", data.table = FALSE)

sp <- read.csv("dat/try/species_list.csv")

trait_ids <- readLines("dat/try/trait_list.txt")

plots <- readRDS("dat/plots_trmm.rds")

ba_mat <- readRDS("dat/ba_clust_mat.rds")

# Remove genera from species list
sp_clean <- sp[unlist(lapply(strsplit(sp$species, " "), length)) > 1,]

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
  filter(species_id %in% sp$id)

# Trait ID match
trait_id_lookup <- try_clean %>%
  dplyr::select(trait_id, trait) %>%
  filter(
    !is.na(trait_id),
    trait_id %in% trait_ids) %>%
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

# Get wood density
wd_sp <- BIOMASS::wdData %>%
  mutate(species = paste(genus, species)) %>%
  filter(
    species %in% sp_clean$species,
    grepl("africa", region, ignore.case = TRUE)) %>%
  group_by(species) %>%
  summarise(val_std = mean(wd, na.rm = TRUE)) %>%
  mutate(trait_id = "wd")

# Bind all traits
traits <- bind_rows(wd_sp, try_species) %>%
  filter(!is.na(val_std)) %>%
  group_by(species, trait_id) %>%
  summarise(val_std = mean(val_std, na.rm = TRUE)) %>%
  spread(trait_id, val_std)

# Create tidy dataframe from ba_mat
ba_df <- ba_mat %>%
  rownames_to_column("plot_cluster") %>%
  gather(species, ba, -plot_cluster) %>%
  filter(ba > 0) %>%
  left_join(., traits, by = "species")

# Calculate community weighted means of traits
plots_cwm_species <- ba_df %>%
  group_by(plot_cluster) %>%
  summarise(
    across(
      all_of(trait_id_lookup$trait_col[trait_id_lookup$trait_col %in% names(.)]),
      ~weighted.mean(.x, ba, na.rm = TRUE),  .names = "{.col}_cwm")) %>%
  as.data.frame(.) %>%
  left_join(., plots, by = "plot_cluster")

# Calculate genus-level CWMs of traits
traits$genus <- unlist(lapply(strsplit(traits$species, " "), "[", 1))

traits_genus <- traits %>%
  dplyr::select(-species) %>%
  group_by(genus) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

ba_df$genus <- unlist(lapply(strsplit(ba_df$species, " "), "[", 1))

ba_df_genus <- ba_df %>%
  dplyr::select(-species) %>%
  group_by(plot_cluster, genus) %>%
  summarise(
    across(
      all_of(trait_id_lookup$trait_col[trait_id_lookup$trait_col %in% names(.)]),
        ~mean(.x, na.rm = TRUE)),
    ba = sum(ba, na.rm = TRUE))

plots_cwm_genus <- ba_df_genus %>%
  group_by(plot_cluster) %>%
  summarise(
    across(
      all_of(trait_id_lookup$trait_col[trait_id_lookup$trait_col %in% names(.)]),
      ~weighted.mean(.x, ba, na.rm = TRUE),  .names = "{.col}_genus_cwm")) %>%
  left_join(., plots_cwm_species, by = "plot_cluster")

# Write to file 
saveRDS(plots_cwm_genus, "dat/plots_try.rds") 

