# Get photoperiod (day length) for each plot across the year
# John Godlee (johngodlee@gmail.com)
# 2020-10-15

# Set working directory

# Packages
library(geosphere)
library(sf)

# Import data
dat <- readRDS("dat/plots.rds")

# get latitude
lat <- st_coordinates(dat)[,2]

# Days of year
doy <- seq(1, 365)

# Tidy dataframe
pp_df <- data.frame(plot_cluster = rep(dat$plot_cluster, each = length(doy)), 
  lat = rep(lat, each = length(doy)),
  doy = doy)

# Calculate photoperiod for each row
pp_df$photoperiod <- daylength(pp_df$lat, pp_df$doy)

# Write file
saveRDS(pp_df, "dat/photoperiod_extract.rds")
