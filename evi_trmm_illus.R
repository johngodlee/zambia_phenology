# Create a plot with a time series of TRMM precipitation data and EVI data
# John L. Godlee (johngodlee@gmail.com)
# Last updated: 2023-09-06

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(RcppRoll)
library(lubridate)

# Import data
trmm_stat <- readRDS("dat/trmm.rds")
trmm_raw <- readRDS("dat/trmm_ts.rds")
modis_stat <- readRDS("dat/modis.rds")

# Sample plots
plot_id <- sample(modis_stat$plot_cluster, 25)

# Filter TRMM data to plot
trmm_fil <- trmm_raw %>% 
  filter(plot_cluster %in% plot_id) 

# Defining rolling period to sum precipitation
roll_period <- 20

# Calculate rolling sum of precipitation
# Create date columns 
trmm_roll <- trmm_fil %>% 
  group_by(plot_cluster) %>% 
  mutate(
    roll = roll_suml(precip, roll_period),
    doy = yday(date),
    doy200 = ifelse(doy > 182, doy - 365, doy),
    year = year(date))

# Filter modis data and make long
modis_fil <- modis_stat %>% 
  mutate(Dormancy2 = MidGreendown + 35) %>% 
  filter(plot_cluster %in% plot_id) %>% 
  pivot_longer(-plot_cluster) %>% 
  filter(name %in% c("Dormancy2", "Greenup", "Maturity", "Senescence"))

# TRMM dates clean
trmm_fil <- trmm_stat %>% 
  filter(plot_cluster %in% plot_id) %>% 
  group_by(plot_cluster) %>% 
  summarise(
    trmm_start = mean(trmm_start, na.rm = TRUE),
    trmm_end = mean(trmm_end, na.rm = TRUE)) %>% 
  pivot_longer(-plot_cluster)

# Plot sample of plots showing TRMM rolling sum MODIS dates, and TRMM dates
pdf(file = "img/modis_trmm_sample.pdf", width = 15, height = 12)
ggplot(trmm_roll, aes(x = doy200, y = roll)) + 
  geom_smooth(aes(group = year), 
    colour = "black", linewidth = 0.5, se = FALSE, span = 0.5) + 
  geom_vline(data = modis_fil, aes(xintercept = value, colour = name)) +
  geom_vline(data = trmm_fil, aes(xintercept = value, colour = name), 
    linetype = 2) +
  facet_wrap(~plot_cluster) + 
  theme_bw() + 
  labs(x = "Days after 1st January", y = "10 day precipitation (mm)")
dev.off()

# Compare methods of estimating start of rainy season
ggplot(trmm_stat, aes(x = trmm_start, y = trmm_start10)) + 
  geom_point() +
  geom_abline(colour = "red", linetype = 2)
# My method produces much early starts of season

# Compare methods of estimating end of rainy season
ggplot(trmm_stat, aes(x = trmm_end, y = trmm_end95)) + 
  geom_point() + 
  geom_abline(colour = "red", linetype = 2)
# My method always produces an end of season later than the q95 method
