# Extract statistics from MODIS 250 m (MOD13Q1) phenology time series 
# John Godlee (johngodlee@gmail.com)
# 2020-08-31

# Packages
library(sf)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(zoo)

source("functions.R")

# Read .csv files per granule 
csv_files <- list.files("dat", "*_extract.csv", full.names = TRUE)
csv_list <- lapply(csv_files, read.csv)

# Check that no plots have been included in both granules 
intersect(unique(csv_list[[1]]$plot_cluster), unique(csv_list[[2]]$plot_cluster))

# Join csv files
evi_ts_df <- do.call(rbind, csv_list)
evi_ts_df$date <- as.Date(evi_ts_df$date)

# Exclude a plot with obviously corrupted data
plot_corrupt <- evi_ts_df %>% 
  filter(date > as.Date("2019-08-01") & 
    date < as.Date("2019-09-01") & 
    evi > 60000000) %>% 
  pull(plot_cluster) %>%
  unique()

evi_ts_clean <- evi_ts_df %>%
  filter(plot_cluster != plot_corrupt)

# Downscale EVI measurements
##' divided by a million
evi_ts_clean$evi <- evi_ts_clean$evi * 0.000001

# Plot all time series
pdf(file = "img/ts.pdf", width = 12, height = 8)
ggplot() +
  geom_path(data = evi_ts_clean, aes(x = date, y = evi, group = plot_cluster),
    alpha = 0.2) +
  scale_x_date(date_labels = "%Y-%b", date_breaks = "3 months") +
  theme_panel() +
  theme(legend.position = "none", 
    axis.text = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    panel.grid.minor = element_blank())
dev.off()

# Decompose annual time series - at September, with 2 month overlap on both ends
seasonGet <- function(x, min_date, max_date, date = "date") {
  # Subset data between max and min dates
  out <- x[x[[date]] > as.Date(min_date) & x[[date]] < as.Date(max_date),]

  # Get the "season" of the data, i.e. the year the measurements start
  out$season <- min(format(out[[date]], "%Y"))

  # Get the days from the min date
  out$days <- as.numeric(out[[date]] - as.Date(min_date))

  # Get year
  out$year <- format(out[[date]], "%Y")

  # Recenter days on start of year in middle of growing season
  out$doy <- as.numeric(out$date - as.Date(paste0(format(as.Date(max_date), "%Y"), "-01-01")))
    

  return(out)
}

evi_ts_list <- split(evi_ts_clean, evi_ts_clean$plot_cluster)

evi_seas_list <- lapply(evi_ts_list, function(x) {
  list(
    s_2015 = seasonGet(x, "2015-07-01", "2016-11-01"),
    s_2016 = seasonGet(x, "2016-07-01", "2017-11-01"),
    s_2017 = seasonGet(x, "2017-07-01", "2018-11-01"),
    s_2018 = seasonGet(x, "2018-07-01", "2019-11-01"),
    s_2019 = seasonGet(x, "2019-07-01", "2020-11-01")
  )
})

evi_clean <- do.call(rbind, 
  do.call(rbind, evi_seas_list)
)

loess_span <- 0.25

write(
  commandOutput(loess_span, "loessSpan"),
  file="out/vars.tex", append=TRUE)

pdf(file = "img/ts_smooth.pdf", width = 12, height = 8)
evi_clean %>%
  filter(plot_cluster %in% 
    unique(.$plot_cluster)[sample(seq(length(unique(.$plot_cluster))), 50)]) %>%
  ggplot(aes(x = doy, y = evi)) + 
    geom_path(aes(group = season), colour = pal[1]) + 
    stat_smooth(method = "loess", span = loess_span, colour = "black") + 
    facet_wrap(~plot_cluster) +
    theme_panel() + 
    theme(legend.position = "none",
      strip.text = element_blank(),
      axis.text.x = element_text(size=12, angle=90, vjust=0.5)) + 
    labs(x = "Days from 1st Jan.", y = "EVI")
dev.off()

pdf(file = "img/ts_season_year.pdf", width = 10, height = 8)
ggplot() + 
  geom_path(data = filter(evi_clean, plot_cluster == "ZIS_1"), 
    aes(x = date, y = evi, colour = season, linetype = year))
dev.off()

# Compute yearly curves for each plot
evi_split <- split(evi_clean, evi_clean$plot_cluster)

pred_data <- seq(min(evi_split[[1]]$doy), max(evi_split[[1]]$doy))

loess_list <- lapply(evi_split, function(x) {
  mod <- loess(evi ~ doy, data = x, span = loess_span)
  pred <- predict(mod, newdata = pred_data) 
  data.frame(pred, doy = pred_data)
})

# Calculate key statistics: 
phen_df <- data.frame(plot_cluster = names(loess_list))

# Max and min
phen_df$max_vi <- unlist(lapply(loess_list, function(x) {max(x[1:355, "pred"])}))
phen_df$min_vi <- unlist(lapply(loess_list, function(x) {min(x[1:355, "pred"])}))

# Create list of minimas by position
min_list <- lapply(seq(length(loess_list)), function(x) {
  # Create zoo ts object
  xz <- as.zoo(loess_list[[x]][["pred"]])

  # Get position of max val
  max_pos <- which(xz == phen_df[x, "max_vi"])

  # Subset to left and right of max val
  left <- xz[1:max_pos]
  right <- xz[max_pos:length(xz)]

  # Take lowest minimum to left and right of max 
  return(c(which(xz == min(left)), which(xz == min(right))))
})

## Subset each loess to within the minimas
loess_fil <- lapply(seq(length(loess_list)), function(x) {
  loess_list[[x]][seq(from = min_list[[x]][1], to = min_list[[x]][length(min_list[[x]])]),]
})

# Average VI (avg_vi)
##' Mean of all values within growing season
phen_df$avg_vi <- unlist(lapply(loess_fil, function(x) {mean(x[["pred"]])}))

# Season start date (s1_start)
##' Date when EVI rises 10% above minimum 

## Extract season start doy
phen_df$s1_start <- unlist(lapply(seq(length(loess_fil)), function(x) {
  min_val <- loess_list[[x]][min_list[[x]][1], "pred"] 
  start_val <- min_val + (min_val * 0.1)
  loess_fil[[x]][min(which(loess_fil[[x]][["pred"]] > start_val)), "doy"]
}))

# Season end date (s1_end)
##' Date after start when EVI falls below 10% above minimum 
phen_df$s1_end <- unlist(lapply(seq(length(loess_fil)), function(x) {
  min_val <- loess_list[[x]][min_list[[x]][1], "pred"] 
  start_val <- min_val + (min_val * 0.1)
  loess_fil[[x]][max(which(loess_fil[[x]][["pred"]] > start_val)), "doy"]
}))

loess_fil <- lapply(seq(length(loess_fil)), function(x) {
  loes <- loess_fil[[x]]
  loes[loes$doy >= phen_df[x, "s1_start"] & loes$doy <= phen_df[x, "s1_end"],]
})

# Season length (s1_length)
##' Days between start and end
phen_df$s1_length <- phen_df$s1_end - phen_df$s1_start

# Greening rate (s1_green_rate)
##' Slope of Linear regression between start (10%) and max
s1_greenup_mod_list <- lapply(seq(length(loess_fil)), function(x) {
  loes <- loess_fil[[x]]
  greenup_fil <- loes[loes[["doy"]] < loes[which(loes[["pred"]] == phen_df[x, "max_vi"]), "doy"],]
  if (nrow(greenup_fil) == 0) {
    return(NA)
  } else {
    return(lm(pred ~ doy, data = greenup_fil))
  }
})

phen_df$s1_green_rate <- unlist(lapply(s1_greenup_mod_list, function(x) {
  if (is.atomic(x)) {
    return(NA)
  } else {
    return(x$coefficients[2])
  }
}))

# Senescence rate (s1_senes_rate)
##' Slope of Linear regression between max and end (10%)
s1_senes_mod_list <- lapply(seq(length(loess_fil)), function(x) {
  loes <- loess_fil[[x]]
  greenup_fil <- loes[loes[["doy"]] > loes[which(loes[["pred"]] == phen_df[x, "max_vi"]), "doy"],]
  if (nrow(greenup_fil) == 0) {
    return(NA)
  } else {
    return(lm(pred ~ doy, data = greenup_fil))
  }
})

phen_df$s1_senes_rate <- unlist(lapply(s1_senes_mod_list, function(x) {
  if (is.atomic(x)) {
    return(NA)
  } else {
    return(x$coefficients[2])
  }
}))

# Cumulative VI (cum_vi)
##' Area under curve between start and end, minus minimum
phen_df$cum_vi <- unlist(lapply(seq(length(loess_fil)), function(x) {
  loess_min <- loess_fil[[x]][["pred"]] - phen_df[x, "min_vi"] 
  sum(diff(loess_fil[[x]][["doy"]])*rollmean(loess_min,2))
}))

# Make a plot which demonstrates all the different numeric statistics
growth_stat_plot <- function(x, raw = FALSE) {
  # Predict values of greenup and senescence rate
  if (class(s1_greenup_mod_list[[x]]) == "lm") {
  greenup_pred <- predict(s1_greenup_mod_list[[x]])
  greenup_pred_df <- data.frame(doy = s1_greenup_mod_list[[x]]$model$doy, 
    pred = greenup_pred)
  }

  if (class(s1_senes_mod_list[[x]]) == "lm") {
  senes_pred <- predict(s1_senes_mod_list[[x]])
  senes_pred_df <- data.frame(doy = s1_senes_mod_list[[x]]$model$doy, 
    pred = senes_pred)
  }

  p <- ggplot()

  if (raw) {
    evi_raw <- evi_clean %>% 
      filter(plot_cluster == names(loess_list[x]))

    p <- p + geom_path(data = evi_raw, aes(x = doy, y = evi, group = season),
      alpha = 0.5) 
  }

  p <- p + 
    geom_path(data = loess_list[[x]], aes(x = doy, y = pred), size = 1.5) + 
    geom_path(data = loess_fil[[x]], aes(x = doy, y = pred), 
      colour = pal[1], size = 1.5)

  if (class(s1_greenup_mod_list[[x]]) == "lm") {
    p <- p + geom_path(data = greenup_pred_df, aes(x = doy, y = pred),
      size = 1.5, colour = pal[2]) 
  }
  if (class(s1_senes_mod_list[[x]]) == "lm") {
    p <- p + geom_path(data = senes_pred_df, aes(x = doy, y = pred), 
      size = 1.5, colour = pal[2]) 
  }
  p + 
    geom_vline(xintercept = pred_data[min_list[[x]]], colour = pal[3] ) + 
    geom_vline(xintercept = loess_list[[x]][which(loess_list[[x]][["pred"]] == phen_df[x, "max_vi"]), "doy"], colour = pal[4]) + 
    theme_panel() + 
    labs(x = "Days from 1st Jan.", y = "EVI") + 
    ylim(10000000, 60000000)
}

ts_stat_plot_list <- lapply(sample(seq(length(loess_list)), 50), growth_stat_plot) 

pdf(file = "img/ts_s1_stats.pdf", width = 20, height = 15)
grid.arrange(grobs = ts_stat_plot_list, ncol = 5)
dev.off()

# Example plot on its own

ts_example_plot <- growth_stat_plot(100, raw = TRUE)

pdf(file = "img/ts_example.pdf", width = 8, height = 6)
ts_example_plot + 
  theme(axis.text = element_text(size = 12), 
    axis.title = element_text(size = 12))
dev.off()

# Load 0.05 degree data to see how they match up
old_dat <- readRDS("dat/plots_vipphen.rds")
old_gather <- old_dat %>%
  dplyr::select(plot_cluster, contains("s1"), vipphen_cum_vi, vipphen_avg_vi) %>%
  st_drop_geometry() %>%
  gather(key, old, -plot_cluster) %>%
  mutate(key = gsub("^vipphen_", "", .$key))

compare <- phen_df %>%
  gather(key, new, -plot_cluster) %>%
  left_join(., old_gather, by = c("plot_cluster", "key")) %>%
  filter(!key %in% c("min_vi", "max_vi"), 
  !(old < 150 & key == "s1_start"))

pdf(file = "img/old_new_compare.pdf", width = 12, height = 8)
ggplot() + 
  geom_point(data = compare, aes(x = old, y = new), 
    alpha = 0.8, colour = "black", fill = pal[5], shape = 21) + 
  geom_smooth(data = compare, aes(x = old, y = new), 
    method = "lm", colour = pal[1]) + 
  facet_wrap(~key, scales = "free") + 
  theme_panel() + 
  labs(x = "MODIS VIPPHEN", y = "MOD13Q1")
dev.off()

phen_all <- old_dat %>%
  left_join(., phen_df, by = "plot_cluster") %>% 
  filter(!is.na(s1_start))

saveRDS(phen_all, "dat/plots_phen.rds")
