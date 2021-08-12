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
library(mgcv)  # gam()
library(gratia)  # derivatives()
library(xtable)
library(RcppRoll)

source("functions.R")

# Import data
plots <- readRDS("dat/plots.rds")
vipphen <- readRDS("dat/vipphen.rds")
evi_ts_df <- readRDS("dat/evi.rds")

# Clean raw ts
evi_ts_clean <- evi_ts_df %>%
  filter(evi > 0) %>%
  filter(plot_cluster %in% plots$plot_cluster)

plots_clean <- plots %>%
  filter(plot_cluster %in% evi_ts_clean$plot_cluster)

# Are all plots in the EVI time series data?
stopifnot(all(plots_clean$plot_cluster %in% evi_ts_clean$plot_cluster))

# Decompose annual time series - at September, with 2 month overlap on both ends
evi_ts_list <- split(evi_ts_clean, evi_ts_clean$plot_cluster)

evi_seas_list <- lapply(evi_ts_list, function(x) {
  list(
    s_2000 = seasonGet(x, "2000-07-01", "2001-11-01"),
    s_2001 = seasonGet(x, "2001-07-01", "2002-11-01"),
    s_2002 = seasonGet(x, "2002-07-01", "2003-11-01"),
    s_2003 = seasonGet(x, "2003-07-01", "2004-11-01"),
    s_2004 = seasonGet(x, "2004-07-01", "2005-11-01"),
    s_2005 = seasonGet(x, "2005-07-01", "2006-11-01"),
    s_2006 = seasonGet(x, "2006-07-01", "2007-11-01"),
    s_2007 = seasonGet(x, "2007-07-01", "2008-11-01"),
    s_2008 = seasonGet(x, "2008-07-01", "2009-11-01"),
    s_2009 = seasonGet(x, "2009-07-01", "2010-11-01"),
    s_2010 = seasonGet(x, "2010-07-01", "2011-11-01"),
    s_2011 = seasonGet(x, "2011-07-01", "2012-11-01"),
    s_2012 = seasonGet(x, "2012-07-01", "2013-11-01"),
    s_2013 = seasonGet(x, "2013-07-01", "2014-11-01"),
    s_2014 = seasonGet(x, "2014-07-01", "2015-11-01"),
    s_2015 = seasonGet(x, "2015-07-01", "2016-11-01"),
    s_2016 = seasonGet(x, "2016-07-01", "2017-11-01"),
    s_2017 = seasonGet(x, "2017-07-01", "2018-11-01"),
    s_2018 = seasonGet(x, "2018-07-01", "2019-11-01"),
    s_2019 = seasonGet(x, "2019-07-01", "2020-11-01")
  )
})

evi_clean <- do.call(rbind, do.call(rbind, evi_seas_list)) %>%
  filter(!is.na(evi))

# Create time series plots
pdf(file = "img/evi_ts.pdf", width = 12, height = 8)
ggplot() +
  geom_path(data = evi_clean, aes(x = date, y = evi, group = plot_cluster),
    alpha = 0.2) +
  scale_x_date(date_labels = "%Y-%b", date_breaks = "3 months") +
  theme_panel() +
  theme(legend.position = "none", 
    axis.text = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    panel.grid.minor = element_blank())
dev.off()

# Compute yearly curves for each plot
evi_split <- split(evi_clean, evi_clean$plot_cluster)

pred_data <- seq(min(evi_split[[1]]$doy), max(evi_split[[1]]$doy))

gam_list <- lapply(evi_split, function(x) {
  series <- x
  mod <- gam(evi ~ s(doy), data = x)
  pred <- predict(mod, newdata = data.frame(doy = pred_data))
  out <- list(series = series, mod = mod, pred = data.frame(pred, doy = pred_data))
})

# Calculate key statistics: 
phen_df <- data.frame(plot_cluster = names(gam_list))

# Max and min EVI
phen_df$max_vi <- unlist(lapply(gam_list, function(x) {max(x[[3]][1:355, "pred"])}))
phen_df$max_vi_date <- unlist(lapply(gam_list, function(x) {
    x[[3]][which.max(x[[3]][1:355, "pred"]), "doy"]
  }))
phen_df$min_vi <- unlist(lapply(gam_list, function(x) {min(x[[3]][1:355, "pred"])}))

# Start of growing season
win_len <- 20

##' Using first derivatives of GAM
##' Filter to before max EVI
##' First rolling period where positive slope of GAM exceeds 50% of max slope
phen_df$s1_start <- unlist(lapply(seq(length(gam_list)), function(x) {
  # Extract mod and predicted values
  mod <- gam_list[[x]][[2]] 
  pred <- gam_list[[x]][[3]]
  
  # First derivative
  der <- as.data.frame(derivatives(mod, type = "forward", 
      newdata = data.frame(doy = pred_data)))

  # Filter to before max VI
  max_date <- pred[pred$pred == phen_df[x, "max_vi"], "doy"]
  der_fil <- der[der$data < max_date, c(3,4)]

  # Find rolling period minimum slopes
  roll <- roll_min(der_fil$derivative, n = win_len, by = 1, align = "left", fill = NA)

  # Find 50% of max slope
  max_slope_half <- max(der_fil$derivative) * 0.5

  # Find first rolling period where all slope exceed 50% of max slope
  sos <- der_fil[roll >= max_slope_half, "data"][1]

  return(sos)
}))

# End of growing season

##' Using first derivatives of GAM
##' Filter to after max EVI
##' Final 1 day period where negative slope of GAM exceeds 50% of GAM
phen_df$s1_end <- unlist(lapply(seq(length(gam_list)), function(x) {
  mod <- gam_list[[x]][[2]] 
  pred <- gam_list[[x]][[3]]

  # First derivative
  der <- as.data.frame(derivatives(mod, type = "backward", 
      newdata = data.frame(doy = pred_data)))

  # Filter to after max VI
  max_date <- pred[pred$pred == phen_df[x, "max_vi"], "doy"]
  der_fil <- der[der$data > max_date, c(3,4)]

  # Find rolling period maximum slopes
  roll <- roll_max(der_fil$derivative, n = win_len, by = 1, align = "right", fill = NA)

  # Find 50% of min slope
  min_slope_half <- min(der_fil$derivative) * 0.5

  # Find last rolling period where all slope below 50% of min slope
  eos <- der_fil[roll <= min_slope_half , "data"]
  eos <- eos[length(eos)]

  return(eos)
}))

# Season length
##' Days between start and end
phen_df$s1_length <- phen_df$s1_end - phen_df$s1_start

# Find bogus values 
phen_df_fil <- phen_df %>%
  filter(s1_start > -200,
    s1_end < 300,
    s1_length < 500,
    max_vi > 1000, 
    min_vi > 500)

# Filter plots to remove bogus values
gam_list_fil <- gam_list[names(gam_list) %in% unique(phen_df_fil$plot_cluster)]

# Subset GAM predicted values to within start and end of season
gam_list_seas <- lapply(seq(length(gam_list_fil)), function(x) {
  pred <- gam_list_fil[[x]][[3]]
  pred <- pred[pred$doy >= phen_df_fil[x, "s1_start"] & 
    pred$doy <= phen_df_fil[x, "s1_end"],]
  gam_list_fil[[x]][["seas"]] <- pred
  return(gam_list_fil[[x]])
})
names(gam_list_seas) <- names(gam_list_fil)

# Average VI 
##' Mean of all values within growing season
phen_df_fil$avg_vi <- unlist(lapply(gam_list_seas, function(x) {
    mean(x[[4]][["pred"]])
}))

# Greening rate (s1_green_rate)
##' Using first derivatives of GAM
##' Filter to before max EVI
##' First rolling period where positive slope of GAM exceeds 50% of max slope
s1_greenup_mod_list <- lapply(seq(length(gam_list_seas)), function(x) {
  mod <- gam_list_seas[[x]][[2]]
  pred <- gam_list_seas[[x]][[3]]

  # Get first derivative of gam
  der <- as.data.frame(derivatives(mod, newdata = data.frame(doy = pred_data)))

  # Filter first derivative to within greening period
  der_fil <- der[der$data >= phen_df_fil[x, "s1_start"] & 
    der$data <= phen_df_fil[x, "max_vi_date"], c(3,4)] 

  # find 50% max slope
  max_slope_half <- max(der_fil$derivative) * 0.5

  # Define start of greening period
  doy_green_start <- round(der_fil[1, "data"])

  # Define end of greening period
  doy_green_end <- first(der_fil[der_fil$derivative <= max_slope_half, "data"])
  
  # Subset gam predicted values
  pred_green <- pred[pred$doy >= doy_green_start & pred$doy <= doy_green_end,]

  # Linear model
  if (all(is.na(pred_green$pred))) {
    mod_green <- NA
  } else {
    mod_green <- lm(pred ~ doy, data = pred_green)
  }

  # Extract slope
  if (is.atomic(mod_green)) {
    slope <- NA
  } else {
    slope <- unname(mod_green$coefficients[2])
  }

  list(mod = mod_green, slope = slope, doy_end = doy_green_end)
})

phen_df_fil$s1_green_rate <- unlist(lapply(s1_greenup_mod_list, `[[`, 2))

# Senescence rate (s1_senes_rate)
##' Second to last -5 to last -5 
s1_senes_mod_list <- lapply(seq(length(gam_list_seas)), function(x) {
  mod <- gam_list_seas[[x]][[2]]
  pred <- gam_list_seas[[x]][[3]]

  # Get first derivative of gam
  der <- as.data.frame(derivatives(mod, newdata = data.frame(doy = pred_data)))

  # Filter first derivative to within senescence period
  der_fil <- der[der$data >= phen_df_fil[x, "max_vi_date"] & 
    der$data < phen_df_fil[x, "s1_end"],] 

  # find 50% min slope
  min_slope_half <- min(der_fil$derivative) * 0.5

  # Define end of senescence period
  doy_senes_end <- der_fil[nrow(der_fil), "data"]

  # Define start of senescence period
  # First time where slope dips below half max negative slope
  doy_senes_start <- der_fil[der_fil$derivative <= min_slope_half, "data"][1]

  # Subset gam predicted values
  pred_senes <- pred[pred$doy >= doy_senes_start & pred$doy <= doy_senes_end,]

  # Linear model
  if (all(is.na(pred_senes$pred))) {
    mod_senes <- NA
  } else {
    mod_senes <- lm(pred ~ doy, data = pred_senes)
  }

  # Extract slope
  if (is.atomic(mod_senes)) {
    slope <- NA
  } else {
    slope <- unname(mod_senes$coefficients[2])
  }
  
  list(mod = mod_senes, slope = slope, doy_start = doy_senes_start)
})

phen_df_fil$s1_senes_rate <- unlist(lapply(s1_senes_mod_list, `[[`, 2))

# Cumulative VI 
##'Area under curve between start and end, minus minimum
phen_df_fil$cum_vi <- unlist(lapply(seq(length(gam_list_seas)), function(x) {
  sum(diff(gam_list_seas[[x]][[4]][["doy"]]) * 
    rollmean(gam_list_seas[[x]][[4]][["pred"]], 2))
}))

phen_all <- plots_clean %>%
  left_join(., phen_df_fil, by = "plot_cluster") %>%
  left_join(., vipphen, by = "plot_cluster") %>%
  filter(vipphen_n_seasons < 2) %>%
  dplyr::select(-vipphen_n_seasons) %>%
  filter(!is.na(vipphen_s1_length)) %>%
  filter(!is.na(s1_start))

# Make a plot which demonstrates all the different numeric phenology stats
stat_plot <- function(x, raw = FALSE, title = TRUE) {
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
    p <- p + geom_path(data = gam_list_seas[[x]][[1]], 
      aes(x = doy, y = evi, group = season), alpha = 0.5) 
  }

  p <- p + 
    geom_path(data = gam_list_seas[[x]][[3]], aes(x = doy, y = pred), size = 1.5) + 
    geom_path(data = gam_list_seas[[x]][[4]], aes(x = doy, y = pred), 
      colour = pal[1], size = 1.5)

  if (class(s1_greenup_mod_list[[x]][[1]]) == "lm") {
    p <- p + 
      geom_path(data = greenup_pred_df, aes(x = doy, y = pred),
        size = 2, lty = "41", linetype = "dotted", colour = pal[2]) +
      geom_vline(xintercept = s1_greenup_mod_list[[x]][[3]],
        colour = "green")
  }
  if (class(s1_senes_mod_list[[x]][[1]]) == "lm") {
    p <- p + 
      geom_path(data = senes_pred_df, aes(x = doy, y = pred), 
        size = 2, lty = "41", linetype = "dotted", colour = pal[2]) +
      geom_vline(xintercept = s1_senes_mod_list[[x]][[3]],
        colour = "green")
  }
  p <- p + 
    geom_vline(xintercept = gam_list_seas[[x]][[3]][
      which(gam_list_seas[[x]][[3]][["pred"]] == phen_df_fil[x, "max_vi"]), "doy"], 
      colour = pal[4]) + 
    geom_vline(xintercept = st_drop_geometry(phen_all[phen_all$plot_cluster == names(gam_list_seas[x]), "s1_start"]) %>% pull(), colour = "blue") + 
    geom_vline(xintercept = st_drop_geometry(phen_all[phen_all$plot_cluster == names(gam_list_seas[x]), "s1_end"]) %>% pull(), colour = "blue") + 
    theme_panel() + 
    labs(x = "Days from 1st Jan.", y = "EVI") 

  if (title) { 
    p <- p + ggtitle(names(gam_list_seas[x]))
  }

  p
}

sam <- sample(seq(length(gam_list_seas)), 50)
ts_stat_plot_list <- lapply(sam, stat_plot) 

# pdf(file = "img/ts_s1_stats.pdf", width = 20, height = 15)
# grid.arrange(grobs = ts_stat_plot_list, ncol = 5)
# dev.off()

# Example plot on its own
ts_example_plot <- stat_plot(598, raw = TRUE, title = FALSE)

pdf(file = "img/ts_example.pdf", width = 8, height = 6)
ts_example_plot + 
  theme(axis.text = element_text(size = 12), 
    axis.title = element_text(size = 12))
dev.off()

# Compare VIPPHEN and 250 m
old_gather <- phen_all %>%
  st_drop_geometry() %>%
  dplyr::select(plot_cluster, starts_with("vipphen"), -vipphen_bg_vi) %>%
  gather(key, old, -plot_cluster) %>%
  mutate(key = gsub("^vipphen_", "", .$key))

compare <- phen_all %>%
  st_drop_geometry() %>%
  dplyr::select(plot_cluster, avg_vi, cum_vi, starts_with("s1")) %>%
  gather(key, new, -plot_cluster) %>%
  left_join(., old_gather, by = c("plot_cluster", "key")) %>%
  filter(!key %in% c("min_vi", "max_vi"), 
  !(old < 150 & key == "s1_start")) %>%
  mutate(key_clean = factor(key, 
      levels = c("avg_vi", "cum_vi", "s1_start", "s1_end", 
        "s1_length", "s1_green_rate", "s1_senes_rate"),
      labels = c("Mean EVI", "Cumulative EVI", "Season start", "Season end",
        "Season length", "Green-up rate", "Senescence rate")))
  
excl <- unique(unname(unlist(list(
  compare[compare$old < 1000 & compare$key == "avg_vi","plot_cluster"],
  compare[compare$old < 600 & compare$key == "cum_vi","plot_cluster"],
  compare[compare$new < 100 & compare$key == "s1_end","plot_cluster"],
  compare[compare$old < 200 & compare$key == "s1_start","plot_cluster"],
  compare[compare$new < 150 & compare$key == "s1_length","plot_cluster"],
  compare[compare$new > 75 & compare$key == "s1_green_rate","plot_cluster"],
  compare[compare$old > 3000 & compare$key == "s1_senes_rate", "plot_cluster"]))))

compare_clean <- compare %>%
  filter(!plot_cluster %in% excl)

compare_list <- split(compare_clean, compare_clean$key_clean)

annot_df <- do.call(rbind, lapply(compare_list, function(x) {
  mod <- lm(new ~ old, data = x)
  summ <- summary(mod)
  pval <- anova(mod)[1,5]
  r2 <- summ$r.squared
  fstat <- summ$fstatistic
  f <- round(fstat[1], 1)
  pfmt <- case_when(
      pval <= 0.05 ~ "p<0.05",
      pval <= 0.01 ~ "p<0.01",
      pval <= 0.001 ~ "p<0.001",
      TRUE ~ paste0("p=", as.character(round(pval, 2))))
  annot <- bquote(R^2~round(r2, 2))
  data.frame(key_clean = unique(x$key_clean), dof = fstat[3], f, r2, pfmt)
  }))

annot_xtable <- xtable(annot_df, 
  label = "annot_df",
  align = "rrcccc",
  display = c("s", "s", "f", "f", "f", "s"),
  digits = c(0, 0, 0, 1, 2, 0),
  caption = "Model fit statistics for comparison of MODIS VIPPHEN and MOD13Q1 products across each of our study sites.")
names(annot_xtable) <- c("Response", "DoF", "F","R\\textsuperscript{2}","Prob.")

fileConn <- file("out/vipphen_compare.tex")
writeLines(print(annot_xtable, include.rownames = FALSE, 
    table.placement = "h",
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)
  
pdf(file = "img/vipphen_compare.pdf", width = 12, height = 8)
ggplot() + 
  geom_point(data = compare_clean, aes(x = old, y = new), 
    alpha = 0.8, colour = "black", fill = pal[5], shape = 21) + 
  geom_smooth(data = compare_clean, aes(x = old, y = new), 
    method = "lm", colour = pal[1]) + 
  facet_wrap(~key_clean, scales = "free") + 
  theme_panel() + 
  labs(x = "MODIS VIPPHEN", y = "MOD13Q1")
dev.off()

# Write data 
saveRDS(phen_all, "dat/plots_phen.rds")
saveRDS(gam_list_seas, "dat/evi_all.rds")

write(
  c(
    commandOutput(win_len, "modisWin"),
    commandOutput(length(excl), "vipphenOutlier")
  ),
  file = "out/modis_extract_vars.tex")
