# Extract TRY database traits - conservative/acquisitive strategies
# John Godlee (johngodlee@gmail.com)
# 2021-02-06

# Community weighted means of functional traits per site, at genus level

# Packages
library(sf)
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)

# Import data
try_dat <- readRDS("dat/traits_sp.rds")

plots <- readRDS("dat/plots_trmm.rds")

# Calculate community weighted means of traits
plots_cwm_species <- try_dat %>%
  group_by(plot_cluster) %>%
  summarise(
    across(
      all_of(trait_id_lookup$trait_col[trait_id_lookup$trait_col %in% names(.)]),
      ~weighted.mean(.x, ba, na.rm = TRUE),  .names = "{.col}_cwm")) %>%
  as.data.frame(.) %>%
  inner_join(., plots, by = "plot_cluster")

# Calculate genus-level CWMs of traits
traits$genus <- unlist(lapply(strsplit(traits$species, " "), "[", 1))

traits_genus <- traits %>%
  dplyr::select(-species) %>%
  group_by(genus) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

try_dat$genus <- unlist(lapply(strsplit(try_dat$species, " "), "[", 1))

try_dat_genus <- try_dat %>%
  dplyr::select(-species) %>%
  group_by(plot_cluster, genus) %>%
  summarise(
    across(
      all_of(trait_id_lookup$trait_col[trait_id_lookup$trait_col %in% names(.)]),
        ~mean(.x, na.rm = TRUE)),
    ba = sum(ba, na.rm = TRUE))

plots_cwm_genus <- try_dat_genus %>%
  group_by(plot_cluster) %>%
  summarise(
    across(
      all_of(trait_id_lookup$trait_col[trait_id_lookup$trait_col %in% names(.)]),
      ~weighted.mean(.x, ba, na.rm = TRUE),  .names = "{.col}_genus_cwm")) %>%
  inner_join(., plots_cwm_species, by = "plot_cluster")

# Write to file 
saveRDS(plots_cwm_genus, "dat/plots_try.rds") 

