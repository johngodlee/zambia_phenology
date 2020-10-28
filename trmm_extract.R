# Extract TRMM rainy data statistics 
# John Godlee (johngodlee@gmail.com)
# 2020-10-13

# Set working directory

# Packages
library(ggplot2)
library(dplyr)
library(sf)
library(mgcv)
library(zoo)
library(gratia)
library(gridExtra)
library(RcppRoll)

source("functions.R")

# Import data
dat <- readRDS("dat/plots_phen.rds")

# Read .csv files per granule 
trmm <- readRDS("dat/trmm.rds")

# Decompose annual time series - at September, with 2 month overlap on both ends
trmm_list <- split(trmm, trmm$plot_cluster)

trmm_seas_list <- lapply(trmm_list, function(x) {
  list(
    s_2015 = seasonGet(x, "2015-07-01", "2016-11-01"),
    s_2016 = seasonGet(x, "2016-07-01", "2017-11-01"),
    s_2017 = seasonGet(x, "2017-07-01", "2018-11-01"),
    s_2018 = seasonGet(x, "2018-07-01", "2019-11-01"),
    s_2019 = seasonGet(x, "2019-07-01", "2020-11-01")
  )
})

trmm_clean <- do.call(rbind, do.call(rbind, trmm_seas_list)) %>%
  filter(!is.na(precip))

saveRDS(trmm_clean, "dat/trmm_ts.rds")

# Plot raw data
pdf(file = "img/trmm_ts.pdf", width = 12, height = 8)
ggplot() +
  geom_path(data = trmm_clean, aes(x = date, y = precip, group = plot_cluster),
    alpha = 0.2) +
  scale_x_date(date_labels = "%Y-%b", date_breaks = "3 months") +
  theme_panel() +
  theme(legend.position = "none", 
    axis.text = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    panel.grid.minor = element_blank())
dev.off()

# Split by plot cluster
trmm_split <- split(trmm_clean, trmm_clean$plot_cluster)

# Generate predict data
pred_data <- seq(min(trmm_split[[1]]$doy), max(trmm_split[[1]]$doy))

# Generate GAMs across maximum precip per 7 days
gam_list <- lapply(trmm_split, function(x) {
  roll <- rollmax(x$precip, k = 7, na.pad = TRUE, align = "right")
  roll[roll == 0] <- NA_real_
  x$roll <- roll
  mod <- gam(roll ~ s(doy), data = x)
  pred <- predict(mod, newdata = data.frame(doy = pred_data))
  pred[pred < 0] <- 0
  out <- list(mod = mod, pred = data.frame(pred, doy = pred_data))
  return(out)
})

# Calculate key statistics: 
trmm_df <- data.frame(plot_cluster = names(gam_list))

# Start of rainy season
win_len <- 20

##' Using first derivatives of GAM
##' Filter to before max EVI
##' First rolling period where positive slope of GAM exceeds 50% of max slope
trmm_df$trmm_start <- unlist(lapply(seq(length(gam_list)), function(x) {
  # Extract mod and predicted values
  mod <- gam_list[[x]][[1]] 
  pred <- gam_list[[x]][[2]] 

  # First derivative
  der <- as.data.frame(derivatives(mod, type = "forward", 
    newdata = data.frame(doy = pred_data)))
  der <- der[,c(3,4)]

  # Find rolling period minimum slopes
  roll <- roll_min(der$derivative, n = win_len, by = 1, align = "left", fill = NA)

  # Find 50% of max slope
  max_slope_half <- max(der$derivative) * 0.5

  # Find first rolling period where all slope exceed 50% of max slope
  sos <- der[roll >= max_slope_half, "data"][1]

  return(sos)
}))

trmm_df$trmm_end <- unlist(lapply(seq(length(gam_list)), function(x) {
  # Extract mod and predicted values
  mod <- gam_list[[x]][[1]] 
  pred <- gam_list[[x]][[2]] 

  # First derivative
  der <- as.data.frame(derivatives(mod, type = "backward", 
    newdata = data.frame(doy = pred_data)))
  der <- der[,c(3,4)]

  # Find rolling period maximum slopes
  roll <- roll_max(der$derivative, n = win_len, by = 1, align = "right", fill = NA)

  # Find 50% of max slope
  min_slope_half <- min(der$derivative) * 0.5

  # Find first rolling period where all slope exceed 50% of max slope
  eos <- der[roll <= min_slope_half, "data"]
  eos <- eos[length(eos)]

  return(eos)
}))

# Season length
##' Days between start and end
trmm_df$trmm_length <- trmm_df$trmm_end - trmm_df$trmm_start

# Cumulative precipitation in wet season
gam_fil <- lapply(seq(length(gam_list)), function(x) {
  mod <- gam_list[[x]][[2]]
  mod[mod$doy >= trmm_df[x, "trmm_start"] & mod$doy <= trmm_df[x, "trmm_end"],]
})

trmm_df$cum_precip_seas <- unlist(lapply(seq(length(gam_fil)), function(x) {
  sum(diff(gam_fil[[x]][["doy"]]) * rollmean(gam_fil[[x]][["pred"]], 2))
}))

# Cumulative precipitation in dry season before wet season (pre-season 90 days)
trmm_df$cum_precip_pre <- unlist(lapply(seq(length(gam_list)), function(x) {
  mod <- gam_list[[x]][[2]]
  season_start <- trmm_df[x, "trmm_start"]
  mod_fil <- mod[mod$doy < season_start & mod$doy >= (season_start - 90),]
  sum(diff(mod_fil[,"doy"]) * rollmean(mod_fil["pred"], 2))
}))

# Cumulative precipitation in dry season before end of season (end season -90 days)
trmm_df$cum_precip_end <- unlist(lapply(seq(length(gam_list)), function(x) {
  mod <- gam_list[[x]][[2]]
  season_end <- trmm_df[x, "trmm_end"]
  mod_fil <- mod[mod$doy < season_end & mod$doy >= (season_end - 90),]
  sum(diff(mod_fil[,"doy"]) * rollmean(mod_fil["pred"], 2))
}))

trmm_df_clean <- dat %>%
  left_join(., trmm_df, by = "plot_cluster") %>%
  filter(!is.na(trmm_start))

gam_list_clean <- gam_list[names(gam_list) %in% unique(trmm_df_clean$plot_cluster)]

# Make a plot which demonstrates all the different numeric phenology stats
stat_plot <- function(x, raw = FALSE) {

  trmm_df_clean <- as.data.frame(st_drop_geometry(trmm_df_clean))
  # Predict values of greenup and senescence rate
  p <- ggplot()

  if (raw) {
    trmm_raw <- trmm_clean %>% 
      filter(plot_cluster == names(gam_list_clean[x]))

    p <- p + geom_path(data = trmm_raw, aes(x = doy, y = precip, group = season),
      alpha = 0.5) 
  }

  p <- p + 
    geom_path(data = gam_list_clean[[x]][[2]], aes(x = doy, y = pred), size = 1.5)

  p + 
    geom_vline(xintercept = trmm_df_clean[
        trmm_df_clean$plot_cluster == names(gam_list_clean[x]), "trmm_start"], 
      colour = "green") + geom_vline(xintercept = trmm_df_clean[
        trmm_df_clean$plot_cluster == names(gam_list_clean[x]), "trmm_end"], 
      colour = "brown") + 
    theme_panel() + 
    labs(x = "Days from 1st Jan.", y = "precip") + 
    ggtitle(names(gam_list_clean[x]))
}

plot_list <- lapply(sample(length(gam_list_clean), 50), stat_plot)

pdf(file = "img/trmm_example.pdf", width = 15, height = 10) 
grid.arrange(grobs = plot_list)
dev.off()

# Write data 
saveRDS(trmm_df_clean, "dat/plots_trmm.rds")

write(
  commandOutput(win_len, "trmmWin"),
  file = "out/trmm_extract_vars.tex")
