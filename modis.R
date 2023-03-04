# Calculate phenological metrics for each site from MODIS MCD12Q2
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-02-24

# Packages
library(dplyr)

source("tex_func.R")

# Read modis data 
out_spread <- readRDS("./dat/modis_raw.rds")
trmm <- readRDS("dat/trmm.rds")

# Calculate length of greening and senescence period
out_spread$green_rate <- out_spread$Maturity - out_spread$Greenup
out_spread$senes_rate <- out_spread$Dormancy - out_spread$Senescence
out_spread$season_length <- out_spread$Dormancy - out_spread$Greenup

# Check all scenese of "best" quality
stopifnot(all(out_spread$QA_Overall == 0 | is.na(out_spread$QA_Overall)))

# Get range of data
modis_range <- range(out_spread$year)

# Calculate average values and lag times
out_avg <- out_spread %>% 
  mutate(
    across(all_of(band_names[1:4]), 
      ~as.Date(.x, origin = "1970-01-01") - as.Date(paste0(year, "-01-01"))),
    across(all_of(band_names[1:4]), 
      ~ifelse(.x > 182, .x-365, .x))) %>% 
  left_join(., trmm, by = c("plot_cluster", "year")) %>%
  mutate(
    start_lag = -(Greenup - trmm_start),
    end_lag = Dormancy - trmm_end) %>% 
  ungroup() %>% 
  dplyr::select(-QA_Overall, -year) %>%
  group_by(plot_cluster) %>% 
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

# Export TeX variables
write(
  c(
    commandOutput(modis_range[1], "modisStart"),
    commandOutput(modis_range[2], "modisEnd")
    ),
  file = "out/modis_vars.tex")


# Write extracted values to file
saveRDS(out_avg, "dat/modis.rds")
