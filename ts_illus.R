# Illustrative plot of pre-rain green-up 
# John Godlee (johngodlee@gmail.com)
# 2021-09-27

# Packages
library(ggplot2)
library(dplyr)

# Import data
trmm_gam <- readRDS("dat/trmm_pred_all.rds")
evi_gam <- readRDS("dat/evi_pred_all.rds")
plot_dat <- readRDS("dat/plots_anal.rds")
plot_trmm_evi_dat <- readRDS("dat/plots_trmm.rds")

# Find a plot which demonstrates pre-rain green-up
plot_dat %>% 
  dplyr::select(plot_cluster, start_lag) %>%
  arrange(desc(start_lag)) %>%
  as.data.frame()

# Make clean dataframe for one plot
chosen_plot <- "ZIS_2170"

trmm_chosen <- trmm_gam[[chosen_plot]]$pred
evi_chosen <- evi_gam[[chosen_plot]]$pred
plot_dat_chosen <- plot_trmm_evi_dat[plot_trmm_evi_dat$plot_cluster == chosen_plot,]

trmm_chosen$key <- "Precipitation (mm)"
evi_chosen$key <- "EVI"

chosen_dat <- bind_rows(trmm_chosen, evi_chosen)

pdf(file = "img/ts_illus.pdf", width = 6, height = 4)
ggplot() +
  geom_line(data = chosen_dat, aes(x = doy, y = pred, colour = key)) + 
  scale_colour_manual(values = c("green", "blue")) + 
  geom_vline(xintercept = plot_dat_chosen$s1_start, colour = "green", linetype = 2) + 
  geom_vline(xintercept = plot_dat_chosen$trmm_start, colour = "blue", linetype = 2) + 
  geom_vline(xintercept = plot_dat_chosen$s1_end, colour = "green", linetype = 2) + 
  geom_vline(xintercept = plot_dat_chosen$trmm_end, colour = "blue", linetype = 2) + 
  facet_wrap(~key, ncol = 1, scales = "free_y") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  labs(x = "DOY", y = "")
dev.off()

