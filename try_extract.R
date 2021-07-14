# Extract TRY database traits - conservative/acquisitive strategies
# John Godlee (johngodlee@gmail.com)
# 2021-02-06

# Community weighted means of functional traits per site, at genus level

# Packages
library(sf)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# Import data
traits_sp <- readRDS("dat/traits_sp.rds")

plots <- readRDS("dat/plots_trmm.rds")

trait_id_lookup <- read.csv("dat/trait_id_lookup.csv")

ba_clust_mat <- readRDS("dat/ba_clust_mat.rds")

# Create genus level means of traits
traits_gen <- traits_sp %>%
  ungroup() %>%
  mutate(genus = gsub(" .*", "", binom)) %>%
  dplyr::select(-binom) %>%
  group_by(genus) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

# Join traits to basal area 
ba_traits <- ba_clust_mat %>% 
  rownames_to_column("plot_cluster") %>%
  gather(binom, ba, -plot_cluster) %>%
  filter(ba > 0) %>%
  mutate(genus = gsub(" .*", "", binom)) %>%
  left_join(., traits_gen, by = "genus")

# What proportion of basal area per site is represented by trait measurements
pdf(file = "img/traits_ecdf.pdf", width = 10, height = 8)
ba_traits %>% 
  group_by(plot_cluster) %>% 
  summarise(across(
      all_of(c("wd", trait_id_lookup$trait_col[trait_id_lookup$trait_col %in% names(.)])),
      ~sum(ba[!is.nan(.x)]) / sum(ba, na.rm = TRUE))) %>%
  gather(trait, prop, -plot_cluster) %>%
  ggplot(.) + 
  geom_line(aes(x = prop, y = 1 - ..y..), 
    stat = "ecdf", pad = FALSE, size = 1) + 
  geom_vline(xintercept = 0.75, linetype = 2, colour = "red") + 
  geom_hline(yintercept = 0.75, linetype = 2, colour = "red") + 
  facet_wrap(~trait) + 
  theme_bw() + 
  labs(x = "Proportion of basal area", y = "Proportion of sites")
dev.off()

# Calculate genus-level CWMs of traits
plots_cwm_genus <- ba_traits %>%
  group_by(plot_cluster) %>%
  summarise(
    across(
      all_of(c("leaf_n_mass", "leaf_p_mass", "sla", "wd")),
      ~weighted.mean(.x, ba, na.rm = TRUE),  .names = "{.col}_genus_cwm")) %>%
  right_join(., plots, by = "plot_cluster")

# Write to file 
saveRDS(plots_cwm_genus, "dat/plots_try.rds") 

