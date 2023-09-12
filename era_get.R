# ERA temperature data processing
# John Godlee (johngodlee@gmail.com)
# 2023-09-12

# Packages
library(dplyr)
library(tidyr)
library(terra)

# Import data
plots <- readRDS("dat/plots.rds")
era_temp <- rast("./dat/ERA_hourly_temp.nc")

# Find daily temperature per plot
era_temp_ext <- terra::extract(era_temp, plots, ID = TRUE)

plot_id_lookup <- data.frame(
  ID = seq_len(nrow(plots)), 
  plot_cluster = plots$plot_cluster)

time_lookup <- data.frame(
  name = as.character(seq_along(unique(time(era_temp)))),
  datetime = unique(time(era_temp)))

era_ext <- era_temp_ext[,c("ID",
    names(era_temp_ext)[grepl("t2m_", names(era_temp_ext))])] %>% 
  pivot_longer(-ID) %>% 
  left_join(., plot_id_lookup, by = "ID") %>%
  mutate(name = gsub("t2m_", "", name)) %>%
  left_join(., time_lookup, by = "name") %>% 
  mutate(date = as.Date(datetime)) %>% 
  group_by(plot_cluster, date) %>% 
  summarise(era_temp = max(value, na.rm = TRUE)) %>% 
  mutate(era_temp = era_temp - 273.15) %>% 
  dplyr::select(plot_cluster, date, era_temp) 

# Write to file
saveRDS(era_ext, "./dat/era_ts.rds")

