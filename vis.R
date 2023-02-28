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

bioclim_zambia <- rast(readRDS("dat/bioclim_zambia.rds"))
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
    diurnal_temp_range_mean = mean(diurnal_temp_range, na.rm = TRUE),
    diurnal_temp_range_sd = sd(diurnal_temp_range, na.rm = TRUE)) %>%
  mutate(
    cluster = as.character(cluster),
    richness = paste0(sprintf("%.0f", richness_median), "(", sprintf("%.0f", richness_iqr), ")"),
    map = paste0(sprintf("%.0f", map_mean), "(", sprintf("%.1f", map_sd), ")"),
    diurnal_temp_range = paste0(sprintf("%.0f", diurnal_temp_range_mean), "(", sprintf("%.1f", diurnal_temp_range_sd), ")")) %>%
  dplyr::select(cluster, n_sites, richness, map, diurnal_temp_range) %>%
  right_join(indval, by = "cluster", multiple = "all") %>%
  mutate(
    species = paste0("\\textit{", species, "}"),
    indval = sprintf("%.3f", indval))

clust_summ[c(rbind(seq(1, nrow(clust_summ), 3), seq(3, nrow(clust_summ), 3))), 
  1:5] <- ""

names(clust_summ) <- c("Cluster", "N sites", "Richness", "MAP", "$\\delta$T", "Species", "Indicator value")

# Export indval table
clust_summ_xtable <- xtable(clust_summ,
  label = "clust_summ",
  align = rep("c", 8),
  display = rep("s", 8),
  caption = "Climatic information and Dufr\\^{e}ne-Legendre indicator species analysis for the vegetation type clusters identified by the PAM algorithm, based on basal area weighted species abundances \\citep{Dufrene1997}. The three species per cluster with the highest indicator values are shown along with other key statistics for each cluster. MAP (Mean Annual Precipitation) and $\\delta$T (Diurnal temperature range) are reported as the mean and 1 standard deviation in parentheses. Species richness is reported as the median and the interquartile range in parentheses.")

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
    cluster = as.character(cluster),
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
phen_bivar_df <- do.call(rbind, lapply(combn(names(resp_lookup), 2, simplify = FALSE),
  function(x) { 
    dat_nogeom <- st_drop_geometry(dat)
    out <- data.frame(dat_nogeom[,x[1]], dat_nogeom[,x[2]], x[1], x[2], 
      dat_nogeom$cluster)
    names(out) <- c("pred", "resp", "x_name", "y_name", "cluster")
    return(out)
  }))

phen_bivar_df$pair_name <- paste(
  phen_bivar_df$x_name,
  phen_bivar_df$y_name,
  sep = "-")

phen_bivar_df$pair_label <- paste(
  "x:", resp_lookup[match(phen_bivar_df$x_name, names(resp_lookup))],
  "\n",
  "y:", resp_lookup[match(phen_bivar_df$y_name, names(resp_lookup))],
  sep = " ")

pdf(file = "img/phen_bivar.pdf", width = 15, height = 8)
ggplot() + 
  geom_point(data = phen_bivar_df, aes(x = pred, y = resp, fill = as.character(cluster)),
    colour = "black", shape = 21) + 
  geom_line(data = phen_bivar_df, aes(x = pred, y = resp),
	stat = "smooth", method = "lm", colour = "black", se = FALSE, linewidth = 1.5) + 
  geom_line(data = phen_bivar_df, aes(x = pred, y = resp, colour = as.character(cluster)), 
	stat = "smooth", method = "lm", se = FALSE, linewidth = 1.5) + 
  facet_wrap(~pair_label, scales = "free", nrow = 3) + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  theme_panel() + 
  labs(x = "", y = "")
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
zambia_trans <- st_transform(zambia, st_crs(greenup_vrt))

# Crop VRTs to Zambia
greenup_crop <- crop(greenup_vrt, zambia_trans)
dormancy_crop <- crop(dormancy_vrt, zambia_trans)

# Find season length
season_length_crop <- dormancy_crop - greenup_crop

# Mask season length layer by Zambia outline
season_length_mask <- mask(season_length_crop, vect(zambia_trans))

# Mask land cover to woodlands and savannas 
lcc_fil <- lcc
lcc_fil[lcc_fil %in% 1:10] <- NA
lcc_wgs84 <- project(lcc_fil, st_crs(zambia)$proj4string)
zambia_sv <- vect(zambia)
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
  geom_sf(data = dat, aes(fill = as.character(cluster)),
    colour = "black", shape = 24, size = 2) +
  theme_panel() + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  labs(x = "", y = ""))
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
    aes(x = mat, y = map, fill = as.character(cluster)),
    shape = 24, size = 2) + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  stat_ellipse(data = dat,
    aes(x = mat, y = map, colour = as.character(cluster)), 
    type = "t", level = 0.95, linewidth = 1.2, show.legend = FALSE) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  theme_panel() + 
  labs(x = expression("MAT" ~ (degree*C)), 
    y = expression("MAP" ~ (mm ~ y^-1))))
dev.off()

# Plot climate and location together
pdf(file = "img/site_map.pdf", width = 15, height = 8)
site_loc + site_clim + 
  plot_layout(guides = "collect")
dev.off()

