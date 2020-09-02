# Testing effect of species composition and richness on greening with Zambia ILUAii data 
# John Godlee (johngodlee@gmail.com)
# 2020-07-29

# Packages
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(sjPlot)
library(spaMM)
library(DHARMa)

source("functions.R")

# Import data 
reli <- readRDS("dat/phen_zam_reliab.rds")
phen_stack <- readRDS("dat/vipphen_stack.rds")

af <- st_read("/Volumes/john/africa_countries/africa.shp")
zambia <- af %>% 
  filter(sov_a3 == "ZMB")

dat <- readRDS("dat/plots_phen.rds")

pcoa_arrows <- readRDS("dat/pcoa_arrows.rds")

# Remove old variables 
dat_clean <- dat %>%
  dplyr::select(-starts_with("vipphen"))

# How many sites are there?
write(
  commandOutput(nrow(dat_clean), "nSites"),
  file="out/vars.tex", append=TRUE)

# Write as shapefile
if (file.exists("dat/shp/loc.shp")) {
  file.remove(list.files("dat/shp", "loc.*", full.names = TRUE))
}
st_write(dat, "dat/shp/loc.shp")

# histogram of raw data
pdf(file =  "img/hist_raw.pdf", width = 12, height = 10)
dat_clean %>% 
  dplyr::select(-plot_id, -plot_id_vec, -plot_cluster) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  gather(variable, value) %>%
  ggplot(aes(x = value)) + 
  geom_histogram(colour = "black", fill = "grey") + 
  facet_wrap(~variable, scales = "free") + 
  labs(x = "", y = "") + 
  theme_bw()
dev.off()

# Create a list of various bivariate relationships
bivar_list <- c(
  "cum_vi ~ richness",
  "cum_vi ~ shannon",
  "cum_vi ~ simpson",
  "cum_vi ~ evenness",
  "cum_vi ~ map",
  "cum_vi ~ mat",

  "s1_length ~ richness",
  "s1_length ~ shannon",
  "s1_length ~ simpson",
  "s1_length ~ evenness",
  "s1_length ~ map",
  "s1_length ~ mat",

  "s1_green_rate ~ richness",
  "s1_green_rate ~ shannon",
  "s1_green_rate ~ simpson",
  "s1_green_rate ~ evenness",
  "s1_green_rate ~ map",
  "s1_green_rate ~ mat",

  "s1_senes_rate ~ richness",
  "s1_senes_rate ~ shannon",
  "s1_senes_rate ~ simpson",
  "s1_senes_rate ~ evenness",
  "s1_senes_rate ~ map",
  "s1_senes_rate ~ mat"
)

# Create models
lm_list <- lapply(bivar_list, function(x) {
  lm(eval(x), data = dat_clean)
})

# Create plots
plot_list <- lapply(bivar_list, function(x) {
  xvar <- unlist(strsplit(x, split = " ~ "))[2]
  yvar <- unlist(strsplit(x, split = " ~ "))[1]
  
  ggplot() + 
    geom_point(data = dat_clean,
	  aes_string(x = xvar, y = yvar),
	  colour = "black", shape = 21, alpha = 0.8) + 
	geom_line(data = dat_clean,
	  aes_string(x = xvar, y = yvar),
	  stat = "smooth", method = "lm", se = FALSE, colour = "blue") + 
	geom_line(data = dat_clean,
	  aes_string(x = xvar, y = yvar), 
	  stat = "smooth", method = "loess", colour = "#DE6400", se = FALSE) + 
	theme_bw()
})

# Arrange on grid
n <- length(plot_list)
n_col <- floor(sqrt(n))

pdf(file =  "img/bivar.pdf", width = 12, height = 10)
do.call("grid.arrange", c(plot_list, ncol = 6))
dev.off()

dat_std <- dat_clean %>% 
  mutate_at(.vars = c(
      "map",
      "mat",
      "diurnal_temp_range",
      "clay",
      "sand",
      "cec",
      "richness",
      "shannon",
      "simpson",
      "evenness", 
      "n_stems_gt5_ha",
      "pcoa_1", "pcoa_2", "pcoa_3", "pcoa_4", "pcoa_5"),
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  dplyr::select(ends_with("_std"), s1_length, cum_vi,
    s1_green_rate, s1_senes_rate, s1_start, geometry) %>%
  st_transform(., UTMProj4("35S")) %>%
  mutate(x = c(unname(st_coordinates(.)[,1])),
    y = c(unname(st_coordinates(.)[,2]))) %>%
  st_drop_geometry()

# Check for collinearity
pdf(file = "img/corrplot.pdf", height = 8, width = 8)
corrPlot(dat_std[,c("richness_std", "evenness_std", "n_stems_gt5_ha_std", 
  "pcoa_1_std", "pcoa_2_std", "pcoa_3_std",
  "map_std", "diurnal_temp_range_std")])
dev.off()
##' None are correlated over r = 0.7, so no serious collinearity

pdf(file = "img/corrplot_response.pdf", height = 8, width = 8)
corrPlot(dat_std[,c("s1_length", "s1_green_rate", "s1_senes_rate", "cum_vi")])
dev.off()

# Linear models ----

# Define model function
phen_mod <- function(var, pre) {
  # Raw data by group
  pdf(file = paste0("img/", pre, "_richness.pdf"), height = 8, width = 10)
  bivar <- ggplot(data = dat_std, 
    aes_string(x = "richness_std", y = var)) + 
    geom_point() + 
    stat_smooth(method = "lm", se = TRUE)
  print(bivar)
  dev.off()

  # Define maximal model
  max_mod <- lm(get(var) ~ richness_std + evenness_std + n_stems_gt5_ha_std + 
    pcoa_1_std + pcoa_2_std + pcoa_3_std + 
    map_std + diurnal_temp_range_std, 
    data = dat_std)

  # Summary
  capture.output(summary(max_mod), file = paste0("out/", pre, "_mod_summary.txt"))

  # QQ plot
  pdf(file = paste0("img/", pre, "_max_mod_qq.pdf"), width = 8, height = 8)
  qqnorm(resid(max_mod))
  qqline(resid(max_mod))
  dev.off()

  # Fixed effect slopes 
  pdf(file = paste0("img/", pre, "_max_mod_slopes.pdf"), width = 12, height = 10)
  fe_slope <- plot_model(max_mod, show.values = TRUE)
  print(fe_slope)
  dev.off()

  # Reduced model comparison 
  max_mod_ml <- lm(get(var) ~ 
    richness_std + evenness_std + n_stems_gt5_ha_std + 
    pcoa_1_std + pcoa_2_std + pcoa_3_std + 
    map_std + diurnal_temp_range_std,
  data = dat_std)

  div_mod_ml <- lm(get(var) ~ 
    richness_std + evenness_std,
    data = dat_std)

  env_mod_ml <- lm(get(var) ~ 
    map_std + mat_std + diurnal_temp_range_std,
    data = dat_std)

  mod_ml_list <- mget(c("max_mod_ml", "div_mod_ml", "env_mod_ml"))

  mod_compare_anova <- eval(parse(text=paste("anova(",
        paste("mod_ml_list[[",1:length(mod_ml_list),"]]",sep="",collapse=","),")")))

  capture.output(mod_compare_anova, file = paste0("out/", pre, "_mod_compare.txt"))
}

# Run function for key responses
phen_mod("s1_length", "l")
phen_mod("s1_green_rate", "r")
phen_mod("s1_start", "s")

# spaMM models with spatial autocorrelation ----

# Reliability mask layer
zambia_clark <- st_transform(zambia, crs(reli))
reli_crop <- crop(reli, zambia_clark)
reli_mask <- mask(reli_crop, zambia_clark)
reli_mask[reli_mask < 2] <- NA
reli_wgs <- projectRaster(reli_mask, phen_wgs_list[[1]])

phen_spamm_mod <- function(max_mod, null_mod, pre) {

  # Model summary
  capture.output(summary(max_mod), 
    file = paste0("out/", pre, "_mod_spamm_summary.txt"))

  pdf(file = paste0("img/", pre, "_max_mod_spamm_slopes.pdf"), width = 6, height = 8)
  spammEffPlot(max_mod)
  dev.off()

  # Estimate degree of spatial autocorrelation
  dd <- dist(dat_std[,c("x","y")])
  mm <- MaternCorr(dd, 
    nu = max_mod$corrPars$`1`$nu, rho = max_mod$corrPars$`1`$rho)
  pdf(file = paste0("img/", pre, "_max_mod_spamm_vario.pdf"), width = 4, height = 4)
  plot(as.numeric(dd), as.numeric(mm), 
    xlab = "Pairwise distance (m)", 
    ylab = "Estimated correlation")
  dev.off()

  sims <- simulateResiduals(max_mod)
  pdf(file = paste0("img/", pre, "_max_mod_spamm_resids.pdf"), width = 8, height = 6)
  plot(sims)
  dev.off()

  # Compare models
  capture.output(anova(max_mod, null_mod), 
    file = paste0("out/", pre, "_mod_spamm_anova.txt"))

  # Map of predicted values
  max_mod_pred <- spammMapExtract(max_mod)
  null_mod_pred <- spammMapExtract(null_mod)
  mod_pred <- cbind(max_mod_pred, null_val = null_mod_pred$val)

  max_mod_raster <- rasterFromXYZ(max_mod_pred, crs = UTMProj4("35S")) 
  zambia_utm <- st_transform(zambia, UTMProj4("35S"))
  max_mod_mask <- mask(max_mod_raster, zambia_utm)

  max_mod_spdf <- as(max_mod_mask, "SpatialPixelsDataFrame")
  max_mod_df <- as.data.frame(max_mod_spdf)
  colnames(max_mod_df) <- c("val", "x", "y")

  pdf(file = paste0("img/", pre, "_mod_spamm_map.pdf"), width = 8, height = 6)
  print(
    ggplot() + 
      stat_contour_filled(data = max_mod_df, aes(x = x, y = y, z = val)) + 
      scale_fill_viridis_d(name = "") + 
      geom_sf(data = st_transform(zambia, UTMProj4("35S")), 
        fill = NA, colour = "black") + 
      geom_point(data = dat_std, aes(x = x, y = y), 
        colour = "black", fill = "white", shape = 24) + 
      theme_classic() + 
      labs(x = "", y = "")
    )
  dev.off()

  if (pre == "l") {
    obs_layer <- phen_wgs_list$`3`
  } 
  else if (pre == "r") {
    obs_layer <- phen_wgs_list$`5`
  } 
  else if (pre == "s") {
    obs_layer <- phen_wgs_list$`1`
  }

  # Compare predicted values of models with observed s1_length across region
  obs_layer_fil <- mask(obs_layer, reli_wgs)
  obs_layer_utm <- projectRaster(obs_layer_fil, crs = UTMProj4("35S"))
  obs_layer_utm_mean <- calc(obs_layer_utm, mean)

  mod_pred$obs_val <- extract(obs_layer_utm_mean, mod_pred[,c(1,2)])

  mod_pred_gather <- gather(mod_pred, model, value, -x, -y, -obs_val) %>%
    filter(obs_val > 0)

  pdf(file = paste0("img/", pre, "_mod_spamm_mod_pred_compare.pdf"), width = 8, height = 8)
  print(
    ggplot() + 
      geom_point(data = mod_pred_gather, aes(x = obs_val, y = value, fill = model),
        colour = "black", shape = 21) + 
      geom_abline(slope = 1, colour = "red", linetype = 2)
    )
  dev.off()

  return(mod_pred)
}

max_mod_spamm_l <- fitme(s1_length ~ richness_std + evenness_std + n_stems_gt5_ha_std + 
    pcoa_1_std + pcoa_2_std + pcoa_3_std + map_std + diurnal_temp_range_std + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_l <- fitme(s1_length ~ map_std + diurnal_temp_range_std + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

max_mod_spamm_r <- fitme(s1_green_rate ~ richness_std + evenness_std + n_stems_gt5_ha_std + 
    pcoa_1_std + pcoa_2_std + pcoa_3_std + map_std + diurnal_temp_range_std + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_r <- fitme(s1_green_rate ~ map_std + diurnal_temp_range_std + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

max_mod_spamm_s <- fitme(s1_start ~ richness_std + evenness_std + n_stems_gt5_ha_std + 
    pcoa_1_std + pcoa_2_std + pcoa_3_std + map_std + diurnal_temp_range_std + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_s <- fitme(s1_start ~ map_std + diurnal_temp_range_std + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

spamm_l_pred <- phen_spamm_mod(max_mod_spamm_l, null_mod_spamm_l, "l")
spamm_r_pred <- phen_spamm_mod(max_mod_spamm_r, null_mod_spamm_r, "r")
spamm_s_pred <- phen_spamm_mod(max_mod_spamm_s, null_mod_spamm_s, "s")

# Plot PCOA with arrows
pdf(file = "img/pcoa_arrows.pdf", width = 10, height = 8)
ggplot() +
  geom_point(data = dat_clean, aes(x = pcoa_1, y = pcoa_2)) +
  geom_segment(data = pcoa_arrows[1:10,],
    aes(xend = Axis.1 / 75, yend = Axis.2 / 75),
    x = 0, y = 0, alpha = 0.7, colour = "red",
    arrow = arrow(length = unit(3, "mm"))) + 
  geom_label_repel(data = pcoa_arrows[1:10,],
    aes(x = Axis.1 / 75, y = Axis.2 / 75, label = species),
    label.padding = unit(0.1, "lines"), size = 3) + 
  scale_colour_viridis(name = "Growth season\nlength (days)") +
  theme_classic() + 
  labs(x = "PCo 1", y = "PCo 2")
dev.off()

# Plot PCOA axes
pcoa_gather <- dat_clean %>% 
  dplyr::select(pcoa_1, pcoa_2, pcoa_3, pcoa_4, pcoa_5, s1_length) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  gather(key, val, -pcoa_1, -s1_length) %>%
  mutate(key = toupper(gsub("oa_", "o ", .$key)))

pdf(file = "img/pcoa.pdf", width = 10, height = 8)
ggplot() + 
  geom_point(data = pcoa_gather, 
    aes(x = pcoa_1, y = val, colour = s1_length)) + 
  scale_colour_viridis(name = "Growth season\nlength (days)") + 
  facet_wrap(~key) + 
  coord_equal() + 
  theme_bw() + 
  labs(x = "PCo 1", y = "")
dev.off()

# Plot location map
s1_length_tile <- as.data.frame(
  as(phen_zam_mean_stack$X3, "SpatialPixelsDataFrame")  # s1_length
)

pdf(file = "img/plot_loc.pdf", height = 10, width = 13)
ggplot() +
  geom_tile(data = s1_length_tile, aes(x = x, y = y, fill = X3)) +
  geom_sf(data = zambia, colour = "black", fill = NA) +
  geom_sf(data = dat_clean, colour = "black", fill = "white", shape = 24) +
  scale_fill_viridis(name = "Growth season\nlength (days)", limits = c(100, 300)) + 
  theme_classic() + 
  labs(x = "", y = "")
dev.off()
