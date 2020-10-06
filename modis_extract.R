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
library(lubridate)
library(mgcv)
library(gratia)

source("functions.R")

# Import data
dat <- readRDS("dat/plots.rds")

# Read .csv files per granule 
csv_files <- list.files("dat", "*_extract.csv", full.names = TRUE)
csv_list <- lapply(csv_files, read.csv)

# Check that no plots have been included in both granules 
any(
  intersect(unique(csv_list[[1]]$plot_cluster), unique(csv_list[[2]]$plot_cluster)),
  intersect(unique(csv_list[[1]]$plot_cluster), unique(csv_list[[3]]$plot_cluster)),
  intersect(unique(csv_list[[1]]$plot_cluster), unique(csv_list[[4]]$plot_cluster)),
  intersect(unique(csv_list[[2]]$plot_cluster), unique(csv_list[[3]]$plot_cluster)),
  intersect(unique(csv_list[[2]]$plot_cluster), unique(csv_list[[4]]$plot_cluster)),
  intersect(unique(csv_list[[3]]$plot_cluster), unique(csv_list[[4]]$plot_cluster))
  )

# Join csv files
evi_ts_df <- do.call(rbind, csv_list)
evi_ts_df$date <- unlist(lapply(evi_ts_df$date, as.Date, 
    tryFormats = c("%Y-%m-%d", "%d/%m/%Y")))

# Decompose annual time series - at September, with 2 month overlap on both ends
seasonGet <- function(x, min_date, max_date, date = "date") {
  # Format date as date
  x[[date]] <- as.Date(x[[date]])

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

evi_ts_list <- split(evi_ts_df, evi_ts_df$plot_cluster)

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
  ) %>%
  filter(!is.na(evi))

saveRDS(evi_clean, "dat/evi_ts.rds")

write(
  commandOutput(loess_span, "loessSpan"),
  file="out/vars.tex", append=TRUE)

# Compute yearly curves for each plot
evi_split <- split(evi_clean, evi_clean$plot_cluster)

pred_data <- seq(min(evi_split[[1]]$doy), max(evi_split[[1]]$doy))

gam_list <- lapply(evi_split, function(x) {
  mod <- gam(evi ~ s(doy), data = x)
  pred <- predict(mod, newdata = data.frame(doy = pred_data))
  out <- list(mod = mod, pred = data.frame(pred, doy = pred_data))
})

# Calculate key statistics: 
phen_df <- data.frame(plot_cluster = names(gam_list))

# Max and min
phen_df$max_vi <- unlist(lapply(gam_list, function(x) {max(x[[2]][1:355, "pred"])}))
phen_df$min_vi <- unlist(lapply(gam_list, function(x) {min(x[[2]][1:355, "pred"])}))


slc <- 5

# Start and end of growing season
##' Using first derivatives of GAM
##' Start = slope first +5
##' End = slope last -5
phen_df$s1_start <- unlist(lapply(gam_list, function(x) {
  mod <- x[[1]] 

  der <- derivatives(mod)

  round(pull(der[which(der$derivative >= slc), "data"])[1])
}))

phen_df$s1_end <- unlist(lapply(gam_list, function(x) {
  mod <- x[[1]] 

  der <- derivatives(mod)

  round(last(pull(der[which(der$derivative <= -slc), "data"])))
}))

# Season length
##' Days between start and end
phen_df$s1_length <- phen_df$s1_end - phen_df$s1_start


# Subset GAM predicted values to within start and end of season
gam_fil <- lapply(seq(length(gam_list)), function(x) {
  mod <- gam_list[[x]][[2]]
  mod[mod$doy >= phen_df[x, "s1_start"] & mod$doy <= phen_df[x, "s1_end"],]
})

# Average VI 
##' Mean of all values within growing season
phen_df$avg_vi <- unlist(lapply(gam_fil, function(x) {mean(x[["pred"]])}))

# Greening rate (s1_green_rate)
##' First +5 to second +5
s1_greenup_mod_list <- lapply(seq(length(gam_list)), function(x) {
  mod <- gam_list[[x]][[1]]
  pred <- gam_list[[x]][[2]]

  # Get first derivative of gam
  der <- as.data.frame(derivatives(mod))

  # Filter first derivative to within season
  der_fil <- der[der$data >= phen_df[x, "s1_start"] & 
    der$data <= phen_df[x, "s1_end"],] 

  # Define start of greening period
  doy_green_start <- round(der_fil[1, "data"])

  # Define end of greening period
  doy_green_end <- first(der_fil[der_fil$derivative <= slc, "data"])
  
  # Subset gam predicted values
  pred_green <- pred[pred$doy >= doy_green_start & pred$doy <= doy_green_end,]

  # Linear model
  mod_green <- lm(pred ~ doy, data = pred_green)

  # Extract slope
  if (is.atomic(mod_green)) {
    slope <- NA
  } else {
    slope <- unname(mod_green$coefficients[2])
  }

  list(mod = mod_green, slope = slope)
})

phen_df$s1_green_rate <- unlist(lapply(s1_greenup_mod_list, `[[`, 2))

# Senescence rate (s1_senes_rate)
##' Second to last -5 to last -5 
s1_senes_mod_list <- lapply(seq(length(gam_list)), function(x) {
  mod <- gam_list[[x]][[1]]
  pred <- gam_list[[x]][[2]]

  # Get first derivative of gam
  der <- as.data.frame(derivatives(mod))

  # Filter first derivative to within season
  der_fil <- der[der$data > phen_df[x, "s1_start"] & 
    der$data < phen_df[x, "s1_end"],] 

  # Define start of senescence period
  doy_senes_start <- last(der_fil[der_fil$derivative >= -slc, "data"])

  # Define end of senescence period
  doy_senes_end <- round(der_fil[nrow(der_fil), "data"])

  # Subset gam predicted values
  pred_senes <- pred[pred$doy >= doy_senes_start & pred$doy <= doy_senes_end,]

  # Linear model
  mod_senes <- lm(pred ~ doy, data = pred_senes)

  # Extract slope
  if (is.atomic(mod_senes)) {
    slope <- NA
  } else {
    slope <- unname(mod_senes$coefficients[2])
  }
  

  list(mod = mod_senes, slope = slope)
})

phen_df$s1_senes_rate <- unlist(lapply(s1_senes_mod_list, `[[`, 2))

# Cumulative VI 
##' Area under curve between start and end, minus minimum
phen_df$cum_vi <- unlist(lapply(seq(length(gam_fil)), function(x) {
  gam_min <- gam_fil[[x]][["pred"]] - phen_df[x, "min_vi"] 
  sum(diff(gam_fil[[x]][["doy"]])*rollmean(gam_min,2))
}))

# Join with original plot data
phen_all <- dat %>%
  left_join(., phen_df, by = "plot_cluster") %>% 
  filter(!is.na(s1_start))

# Write data 
saveRDS(phen_all, "dat/plots_phen.rds")

# Make a plot which demonstrates all the different numeric phenology stats
growth_stat_plot <- function(x, raw = FALSE) {
  # Predict values of greenup and senescence rate
  if (class(s1_greenup_mod_list[[x]][[1]]) == "lm") {
  greenup_pred <- predict(s1_greenup_mod_list[[x]][[1]])
  greenup_pred_df <- data.frame(doy = s1_greenup_mod_list[[x]][[1]]$model$doy, 
    pred = greenup_pred)
  }

  if (class(s1_senes_mod_list[[x]][[1]]) == "lm") {
  senes_pred <- predict(s1_senes_mod_list[[x]][[1]])
  senes_pred_df <- data.frame(doy = s1_senes_mod_list[[x]][[1]]$model$doy, 
    pred = senes_pred)
  }

  p <- ggplot()

  if (raw) {
    evi_raw <- evi_clean %>% 
      filter(plot_cluster == names(gam_list[x]))

    p <- p + geom_path(data = evi_raw, aes(x = doy, y = evi, group = season),
      alpha = 0.5) 
  }

  p <- p + 
    geom_path(data = gam_list[[x]][[2]], aes(x = doy, y = pred), size = 1.5) + 
    geom_path(data = gam_fil[[x]], aes(x = doy, y = pred), 
      colour = pal[1], size = 1.5)

  if (class(s1_greenup_mod_list[[x]][[1]]) == "lm") {
    p <- p + geom_path(data = greenup_pred_df, aes(x = doy, y = pred),
      size = 2, linetype = "dotdash", colour = pal[2]) 
  }
  if (class(s1_senes_mod_list[[x]][[1]]) == "lm") {
    p <- p + geom_path(data = senes_pred_df, aes(x = doy, y = pred), 
      size = 2, linetype = "dotdash", colour = pal[2]) 
  }
  p + 
    geom_vline(xintercept = gam_list[[x]][[2]][which(gam_list[[x]][[2]][["pred"]] == phen_df[x, "max_vi"]), "doy"], colour = pal[4]) + 
    theme_panel() + 
    labs(x = "Days from 1st Jan.", y = "EVI") + 
    ggtitle(names(gam_list[x]))
}

#sam <- sample(seq(length(gam_list)), 50)
ts_stat_plot_list <- lapply(sam, growth_stat_plot) 

pdf(file = "img/ts_s1_stats.pdf", width = 20, height = 15)
grid.arrange(grobs = ts_stat_plot_list, ncol = 5)
dev.off()

# Example plot on its own
ts_example_plot <- growth_stat_plot(598, raw = TRUE)

pdf(file = "img/ts_example.pdf", width = 8, height = 6)
ts_example_plot + 
  theme(axis.text = element_text(size = 12), 
    axis.title = element_text(size = 12))
dev.off()



