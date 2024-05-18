# Prepare an anonymous dataset with enough info for analyses
# John L. Godlee (johngodlee@gmail.com)
# Last updated: 2024-05-08

# Packages
library(dplyr)
library(sf)

# Import data
plots <- readRDS("dat/plots.rds")
div <- readRDS("dat/div.rds")
stat_all <- readRDS("dat/stat_all.rds")
trees <- readRDS("dat/trees.rds")

# Clean data used for modelling
stat_clean <- plots %>% 
  st_drop_geometry() %>% 
  inner_join(., div, by = "plot_cluster") %>% 
  inner_join(., stat_all, by = "plot_cluster") %>%
  dplyr::select(
    plot_cluster,
    plot_id,
    year,
    EVI_Area,
    season_length, 
    start_lag,
    end_lag,
    green_rate,
    senes_rate,
    cum_precip_end,
    cum_precip_pre,
    cum_precip_seas,
    cum_temp_end,
    cum_temp_pre,
    cum_temp_seas,
    eff_rich, 
    diam_quad_mean, 
    Detarioideae)

# Write to file
write.csv(stat_clean, "./dat/dat_anon/dat_anon.csv", row.names = FALSE)
