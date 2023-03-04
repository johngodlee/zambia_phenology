# Visualisations and descriptive tables of data
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-02-25

# Packages
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(xtable)
library(shades)
library(tidyterra)
library(ggnewscale)
library(scico)

source("./plot_func.R")

# Import data 
plots <- readRDS("./dat/plots.rds")
div <- readRDS("./dat/div.rds")
modis <- readRDS("./dat/modis.rds")
bioclim <- readRDS("./dat/bioclim.rds")
indval <- readRDS("./dat/indval.rds")

bioclim_zambia <- readRDS("dat/bioclim_zambia.rds")
af <- st_read("dat/africa_countries/africa.shp")
lcc <- rast("dat/zambia_landcover/lcc.tif")

# Combine plots dataframes
dat <- plots %>% 
  inner_join(., div, by = "plot_cluster") %>% 
  inner_join(., modis, by = "plot_cluster") %>%
  inner_join(., bioclim, by = "plot_cluster")

# Check only plots used for analysis in plots dataframe
stopifnot(nrow(dat) == nrow(div))

# Create density distributions of season start and end dates
start_dens_plot <- dat %>%
  dplyr::select(Greenup, trmm_start, cluster) %>%
  st_drop_geometry() %>% 
  pivot_longer(
    -cluster,
    names_to = "variable",
    values_to = "value") %>%
  ggplot(., aes(x = value, colour = variable)) + 
    geom_density() + 
    scale_colour_manual(name = "", 
      labels = c("Growth season start", "Rainy season start"), 
      values = pal[c(2,1)]) + 
    facet_wrap(~cluster) + 
    theme_panel() + 
    theme(legend.position = "bottom") + 
    labs(x = "DOY", y = "Frequency")

end_dens_plot <- dat %>%
  dplyr::select(Dormancy, trmm_end, cluster) %>%
  st_drop_geometry() %>% 
  pivot_longer(
    -cluster,
    names_to = "variable",
    values_to = "value") %>%
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
start_dens_plot + end_dens_plot + plot_layout(ncol = 1)
dev.off()

dat$cluster <- factor(dat$cluster,
      levels = names(clust_lookup),
      labels = clust_lookup)

# Create table of species indicators and climatic data per cluster
clust_summ <- dat %>% 
  st_drop_geometry() %>% 
  group_by(cluster) %>%
  summarise(
    richness_median = median(richness, na.rm = TRUE),
    richness_iqr = (quantile(richness, 0.75) - quantile(richness, 0.25)),
    n_sites = as.character(n()),
    map_mean = mean(map, na.rm = TRUE),
    map_sd = sd(map, na.rm = TRUE),
    mat_mean = mean(mat, na.rm = TRUE),
    mat_sd = sd(mat, na.rm = TRUE)) %>%
  mutate(
    richness = paste0(sprintf("%.0f", richness_median), "(", sprintf("%.0f", richness_iqr), ")"),
    map = paste0(sprintf("%.0f", map_mean), "(", sprintf("%.1f", map_sd), ")"),
    mat = paste0(sprintf("%.0f", mat_mean), "(", sprintf("%.1f", mat_sd), ")")) %>%
  dplyr::select(cluster, n_sites, richness, map, mat) %>%
  right_join(indval, by = "cluster", multiple = "all") %>%
  mutate(
    species = paste0("\\textit{", species, "}"),
    indval = sprintf("%.3f", indval),
    cluster = paste("{\\multirow{3}{*}{\\makecell[c]{",
      gsub("\\s", "\\\\\\\\", cluster),
      "}}}"))

clust_summ[c(rbind(seq(2, nrow(clust_summ), 3), seq(3, nrow(clust_summ), 3))), 
  1] <- ""
clust_summ[c(rbind(seq(1, nrow(clust_summ), 3), seq(3, nrow(clust_summ), 3))), 
  2:5] <- ""

names(clust_summ) <- c("Cluster", "N sites", "Richness", "MAP", "MAT", "Indicator species", "Indicator value")

# Export indval table
clust_summ_xtable <- xtable(clust_summ,
  label = "clust_summ",
  align = rep("c", 8),
  display = rep("s", 8),
  caption = "Climatic information and Dufr\\^{e}ne-Legendre indicator species analysis for the vegetation type clusters identified by the PAM algorithm, based on basal area weighted species abundances \\citep{Dufrene1997}. The three species per cluster with the highest indicator values are shown along with other key statistics for each cluster. MAP (Mean Annual Precipitation) and MAT (Mean Annual Temperature) are reported as the mean and 1 standard deviation in parentheses. Species richness is reported as the median and the interquartile range in parentheses.")

fileConn <- file("out/clust_summ.tex")
writeLines(print(clust_summ_xtable, include.rownames = FALSE,
    table.placement = "H",
    hline.after = c(-1,0,seq(from = 3, by = 3, length.out = length(unique(indval$cluster)))),
    sanitize.text.function = function(x) {x}),
  fileConn)
close(fileConn)

# Density plots of phenological metrics per cluster
pdf(file = "img/phen_dens_clust.pdf", width = 12, height = 10)
dat %>%
  st_drop_geometry() %>% 
  dplyr::select(names(resp_lookup), plot_cluster, cluster) %>%
  pivot_longer(
    -c(cluster, plot_cluster),
    names_to = "variable",
    values_to = "value") %>%
  mutate(
    cluster = cluster,
    variable = factor(variable, 
      levels = names(resp_lookup)[c(1,3,5,2,4,6)],
      labels = resp_lookup[c(1,3,5,2,4,6)])) %>%
  ggplot(., aes(x = value, group = cluster, colour = cluster)) + 
  geom_density(linewidth = 1.5) + 
  facet_wrap(~variable, scales = "free") + 
  labs(x = "", y = "") +
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  theme_panel()
dev.off()
    
# Scatter plots comparing each phenological metric
bivar_plot_list <- apply(combn(names(resp_lookup), 2), 2, function(x) {
  ggplot(data = st_drop_geometry(dat), 
      aes(x = .data[[x[1]]], y = .data[[x[2]]])) + 
    geom_point(aes(fill = cluster), 
      colour = "black", shape = 21) + 
    geom_line(aes(colour = cluster),
      stat = "smooth", method = "lm", se = FALSE, linewidth = 1.5) + 
    geom_line(stat = "smooth", method = "lm", se = FALSE, linewidth = 1.5) + 
    scale_fill_manual(name = "Cluster", values = clust_pal) + 
    scale_colour_manual(name = "Cluster", values = clust_pal) + 
    theme_panel() + 
    labs(
      x = resp_plot_axes[names(resp_plot_axes) == x[1]], 
      y = resp_plot_axes[names(resp_plot_axes) == x[2]])
})

pdf(file = "img/phen_bivar.pdf", width = 10, height = 12)
wrap_plots(bivar_plot_list) +  
  plot_layout(
    ncol = 3, 
    guides = "collect")  & theme(legend.position = 'bottom')
dev.off()

# Get MODIS HDF file names for 2014 scenes
modis_files <- list.files("./dat/MCD12Q2_2014", pattern = "*.hdf", full.names = TRUE)

# Get band names for Greenup and Dormancy
greenup_band_list <- paste0("HDF4_EOS:EOS_GRID:", modis_files, 
  ":MCD12Q2:", "Greenup")
dormancy_band_list <- paste0("HDF4_EOS:EOS_GRID:", modis_files, 
  ":MCD12Q2:", "Dormancy")

# Create VRTs for Greenup and Dormancy
greenup_vrt <- vrt(greenup_band_list)[[1]]
dormancy_vrt <- vrt(dormancy_band_list)[[1]]

# Project Zambia outline to same CRS as VRTs
zambia <- af[af$sov_a3 == "ZMB",]
zambia_trans <- st_transform(zambia, st_crs(greenup_vrt))

# Crop VRTs to Zambia
# greenup_crop <- terra::crop(greenup_vrt, zambia_trans)
# dormancy_crop <- terra::crop(dormancy_vrt, st_geometry(zambia_trans))

# saveRDS(greenup_crop, "./dat/MCD12Q2_2014/greenup_crop.rds")
# saveRDS(dormancy_crop, "./dat/MCD12Q2_2014/dormancy_crop.rds")

greenup_crop <- readRDS("./dat/MCD12Q2_2014/greenup_crop.rds")
dormancy_crop <- readRDS("./dat/MCD12Q2_2014/dormancy_crop.rds")

# Find season length
season_length_crop <- dormancy_crop - greenup_crop

# Mask season length layer by Zambia outline
season_length_mask <- mask(season_length_crop, vect(zambia_trans))

# Mask land cover to woodlands and savannas 
lcc_fil <- lcc
lcc_fil[lcc_fil %in% 1:10] <- NA
lcc_wgs84 <- project(lcc_fil, crs(zambia_trans))
zambia_sv <- vect(zambia_trans)
lcc_mask <- mask(crop(lcc_wgs84, zambia_sv), zambia_sv)

# Aggregate to same resolution as phenology
lcc_agg <- resample(lcc_mask, season_length_mask)

# Mask phenology by land cover
season_length_mask_veg <- mask(season_length_mask, lcc_agg, inverse = TRUE)

# Plot location map
pdf(file = "img/site_loc.pdf", height = 8, width = 10)
(site_loc <- ggplot() +
  geom_spatraster(data = season_length_mask_veg) + 
  scale_fill_gradient(name = "Season length\n(days)", 
    low = "black" , high = pal[6], limits = c(100, 300), na.value = NA) + 
  new_scale_fill() +
  geom_sf(data = zambia, colour = "black", fill = NA) +
  geom_sf(data = dat, aes(fill = cluster),
    colour = "black", shape = 24, size = 2) +
  theme_panel() + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  labs(x = "", y = "") + 
  ggtitle("Geographic space"))
dev.off()

# Plots in climate space
bioclim_val <- as.data.frame(values(bioclim_zambia)) %>% 
  filter(!is.na(wc2.1_30s_bio_1))

pdf(file = "img/site_clim.pdf", width = 10, height = 8)
(site_clim <- ggplot() + 
  geom_bin2d(data = bioclim_val, 
    aes(x = wc2.1_30s_bio_1, y = wc2.1_30s_bio_12, fill = after_stat(count)), 
    bins = 100) +
  scale_fill_scico(name = "Pixel density", palette = "bamako",
     trans = "log", breaks = c(1, 10, 100, 1000, 10000)) + 
  new_scale_fill() + 
  geom_point(data = dat, 
    aes(x = mat, y = map, fill = cluster),
    shape = 24, size = 2) + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  stat_ellipse(data = dat,
    aes(x = mat, y = map, colour = cluster), 
    type = "t", level = 0.95, linewidth = 1.2, show.legend = FALSE) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  theme_panel() + 
  labs(x = expression("MAT" ~ (degree*C)), 
    y = expression("MAP" ~ (mm ~ y^-1))) + 
  ggtitle("Climate space"))
dev.off()

# Plot climate and location together
pdf(file = "img/site_map.pdf", width = 15, height = 8)
site_loc + site_clim + 
  plot_layout(guides = "collect")
dev.off()

