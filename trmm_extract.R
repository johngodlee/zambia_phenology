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

source("functions.R")

# Import data
dat <- readRDS("dat/plots_vipphen.rds")

# Read .csv files per granule 
trmm <- read.csv("dat/trmm_extract.csv")

# Decompose annual time series - at September, with 2 month overlap on both ends
trmm_list <- split(trmm, trmm$plot_cluster)

trmm_seas_list <- lapply(trmm_list, function(x) {
  list(
    s_2015 = seasonGet(x, "2015-12-31", "2016-12-31"),
    s_2016 = seasonGet(x, "2016-12-31", "2017-12-31"),
    s_2017 = seasonGet(x, "2017-12-31", "2018-12-31"),
    s_2018 = seasonGet(x, "2018-12-31", "2019-12-31"),
    s_2019 = seasonGet(x, "2019-12-31", "2020-12-31")
  )
})

trmm_clean <- do.call(rbind, do.call(rbind, trmm_seas_list)) %>%
  filter(!is.na(precip))

saveRDS(trmm_clean, "dat/trmm_ts.rds")

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

trmm_split <- split(trmm_clean, trmm_clean$plot_cluster)

pred_data <- seq(min(trmm_split[[1]]$doy), max(trmm_split[[1]]$doy))

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

# Define change in slope parameter
slc <- 0.05

# Start and end of growing season
##' Using first derivatives of GAM
##' Start = slope first +5
##' End = slope last -5
trmm_df$trmm_start <- unlist(lapply(gam_list, function(x) {
  mod <- x[[1]] 

  der <- derivatives(mod, type = "forward")

  der_fil <- der[der$data > 100,]

  round(pull(der_fil[which(der_fil$derivative >= slc), "data"])[1])
}))

trmm_df$trmm_end <- unlist(lapply(gam_list, function(x) {
  mod <- x[[1]] 

  der <- derivatives(mod, type = "backward")

  der_fil <- der[der$data < 350,]

  round(last(pull(der_fil[which(der_fil$derivative <= -slc), "data"])))
}))

# Season length
##' Days between start and end
trmm_df$trmm_length <- trmm_df$trmm_end - trmm_df$trmm_start

trmm_df_clean <- trmm_df %>%
  left_join(., dat, by = "plot_cluster") %>%
  filter(!is.na(trmm_start))

# Write data 
saveRDS(trmm_df_clean, "dat/plots_trmm.rds")

