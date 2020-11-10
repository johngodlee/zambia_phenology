# Pre-processing analysis
# John Godlee (johngodlee@gmail.com)
# 2020-07-29

# Packages
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggnewscale)
library(shades)
library(gridExtra)
library(raster)

source("functions.R")

# Import data 
dat <- readRDS("dat/plots_div.rds")

dat <- st_as_sf(dat)

phen_stack <- readRDS("dat/vipphen_stack.rds")

af <- st_read("dat/africa_countries/africa.shp")
zambia <- af %>% 
  filter(sov_a3 == "ZMB")

# Calculate some statistics
dat_clean <- dat %>%
  mutate(cluster = factor(cluster, 
      labels = clust_lookup[1:length(unique(dat$cluster))]),
    start_lag = -(s1_start - trmm_start),
    end_lag = s1_end - trmm_end)
  
start_dens_plot <- dat_clean %>%
  dplyr::select(s1_start, trmm_start, cluster) %>%
  st_drop_geometry() %>%
  gather(variable, value, -cluster) %>%
  ggplot(., aes(x = value, colour = variable)) + 
    geom_density() + 
    scale_colour_manual(name = "", 
      labels = c("Growth season start", "Rainy season start"), 
      values = pal[c(2,1)]) + 
    facet_wrap(~cluster) + 
    theme_panel() + 
    theme(legend.position = "bottom") + 
    labs(x = "DOY", y = "Frequency")

end_dens_plot <- dat_clean %>%
  dplyr::select(s1_end, trmm_end, cluster) %>%
  st_drop_geometry() %>%
  gather(variable, value, -cluster) %>%
  ggplot(., aes(x = value, colour = variable)) + 
    geom_density() + 
    scale_colour_manual(name = "", 
      labels = c("Growth season end", "Rainy season end"), 
      values = pal[c(2,1)]) + 
    facet_wrap(~cluster) + 
    theme_panel() + 
    theme(legend.position = "bottom") + 
    labs(x = "DOY", y = "Frequency")

pdf(file = "img/dens_lag.pdf", width = 10, height = 12)
grid.arrange(grobs = list(start_dens_plot, end_dens_plot), ncol = 1)
dev.off()

# How many sites are there?
write(
  commandOutput(nrow(dat_clean), "nSites"),
  file = "out/analysis_vars.tex")

# Density plots of phenological metrics per cluster
pdf(file =  "img/phen_dens_clust.pdf", width = 12, height = 10)
dat_clean %>%
  dplyr::select(names(resp_lookup), cluster) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  gather(variable, value, -cluster) %>%
  dplyr::select(variable, value, cluster) %>% 
  mutate(variable = factor(variable, 
      levels = names(resp_lookup)[c(1,3,5,2,4,6)],
      labels = resp_lookup[c(1,3,5,2,4,6)])) %>%
  ggplot(., aes(x = value, colour = cluster)) + 
  geom_density(size = 1.5) + 
  facet_wrap(~variable, scales = "free") + 
  labs(x = "", y = "") +
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  theme_panel()
dev.off()

# Create bivariate relationships plot
bivar_list <- c(
  "cum_vi ~ eff_rich",
  "cum_vi ~ n_stems_gt10_ha",
  "cum_vi ~ map",
  "cum_vi ~ diurnal_temp_range",

  "s1_length ~ eff_rich",
  "s1_length ~ n_stems_gt10_ha",
  "s1_length ~ map",
  "s1_length ~ diurnal_temp_range",

  "s1_green_rate ~ eff_rich",
  "s1_green_rate ~ n_stems_gt10_ha",
  "s1_green_rate ~ map",
  "s1_green_rate ~ diurnal_temp_range",

  "s1_senes_rate ~ eff_rich",
  "s1_senes_rate ~ n_stems_gt10_ha",
  "s1_senes_rate ~ map",
  "s1_senes_rate ~ diurnal_temp_range",

  "start_lag ~ eff_rich",
  "start_lag ~ n_stems_gt10_ha",
  "start_lag ~ map",
  "start_lag ~ diurnal_temp_range", 

  "end_lag ~ eff_rich",
  "end_lag ~ n_stems_gt10_ha",
  "end_lag ~ map",
  "end_lag ~ diurnal_temp_range"
)

bivar_df <- as.data.frame(do.call(rbind, lapply(bivar_list, function(x) {
  x_var <- sym(unlist(strsplit(x, split = " ~ "))[2])
  y_var <- sym(unlist(strsplit(x, split = " ~ "))[1])

  dat_clean %>% 
    dplyr::select(!!x_var, !!y_var, cluster) %>%
    st_drop_geometry() %>%
    rename(pred = !!x_var, resp = !!y_var) %>%
    mutate(x = as.character(x_var), 
      y = as.character(y_var))
})))

bivar_df$x <- factor(bivar_df$x, levels = names(pred_lookup))
bivar_df$y <- factor(bivar_df$y, levels = names(resp_lookup))

pdf(file = "img/bivar.pdf", width = 15, height = 10)
ggplot() + 
  geom_point(data = bivar_df, aes(x = pred, y = resp, fill = cluster), 
	colour = "black", shape = 21) +
  geom_line(data = bivar_df, aes(x = pred, y = resp),
	stat = "smooth", method = "lm", colour = "black", se = FALSE, size = 1.5) + 
  geom_line(data = bivar_df, aes(x = pred, y = resp, colour = cluster), 
	stat = "smooth", method = "lm", se = FALSE) + 
  facet_grid(y~x, scales = "free", 
    labeller = labeller(y = resp_lookup, x = pred_lookup)) +  
  scale_fill_manual(name = "", values = clust_pal) + 
  scale_colour_manual(name = "", values = brightness(clust_pal, 0.5)) + 
  theme_panel() + 
  labs(x = "", y = "")
dev.off()

# Standardise variables
dat_std <- dat_clean %>% 
  mutate_at(.vars = names(pred_lookup)[which(names(pred_lookup) != "cluster")],
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  dplyr::select(ends_with("_std"), cluster, names(resp_lookup), geometry) %>%
  rename_at(.vars = vars(ends_with("_std")), 
    .funs = list(~gsub("_std", "", .))) %>%
  st_transform(., UTMProj4("35S")) %>%
  mutate(x = c(unname(st_coordinates(.)[,1])),
    y = c(unname(st_coordinates(.)[,2]))) %>%
  st_drop_geometry() %>%
  drop_na()

# Check for collinearity
pdf(file = "img/corrplot.pdf", height = 8, width = 8)
corrPlot(dat_std[,names(pred_lookup)[which(names(pred_lookup) != "cluster")]]) + 
  scale_x_discrete(labels = pred_lookup) + 
  scale_y_discrete(labels = pred_lookup)
dev.off()
##' None are correlated over r = 0.7, so no serious collinearity

# Plot location map
s1_length_tile <- as.data.frame(
  as(phen_stack$X3, "SpatialPixelsDataFrame")  # s1_length
)

pdf(file = "img/plot_loc.pdf", height = 8, width = 10)
ggplot() +
  geom_tile(data = s1_length_tile, aes(x = x, y = y, fill = X3)) +
  scale_fill_gradient(name = "Season length\n(days)", low = "black" , high = pal[6], 
    limits = c(100, 300)) + 
  new_scale_fill() +
  geom_sf(data = zambia, colour = "black", fill = NA) +
  geom_sf(data = dat_clean, aes(fill = cluster), 
    colour = "black", shape = 24, size = 3) +
  theme_panel() + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  labs(x = "", y = "")
dev.off()

# Save data ready for models
saveRDS(dat_std, "dat/plots_anal.rds")

