# Testing effect of species composition and richness on greening with Zambia ILUAii data 
# John Godlee (johngodlee@gmail.com)
# 2020-07-29

# Packages
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(gridExtra)
library(sjPlot)
library(spaMM)
library(DHARMa)
library(hier.part)
library(xtable)
library(raster)

source("functions.R")

# Import data 
phen_stack <- readRDS("dat/vipphen_stack.rds")

af <- st_read("/Volumes/john/africa_countries/africa.shp")
zambia <- af %>% 
  filter(sov_a3 == "ZMB")

dat <- readRDS("dat/plots_phen.rds")

pcoa_arrows <- readRDS("dat/pcoa_arrows.rds")

# Define variable name translation lookup
pred_lookup <- c("Richness", "Evenness", "Stem density",
    "PCoA 1", "PCoA 2", "PCoA 3", "MAP", "Diurnal dT")
names(pred_lookup) <- c("richness", "evenness", "n_stems_gt5_ha", 
  "pcoa_1", "pcoa_2", "pcoa_3", "map", "diurnal_temp_range")

resp_lookup <- c("Cumulative EVI", "Season length", 
  "Greening rate", "Senescence rate", "Season start")
names(resp_lookup) <- c("cum_vi", "s1_length", 
  "s1_green_rate", "s1_senes_rate", "s1_start")

# Remove old variables 
dat_clean <- dat %>%
  dplyr::select(-starts_with("vipphen"))

# How many sites are there?
write(
  commandOutput(nrow(dat_clean), "nSites"),
  file="out/vars.tex", append=TRUE)

# Write as shapefile
#if (file.exists("dat/shp/loc.shp")) {
#  file.remove(list.files("dat/shp", "loc.*", full.names = TRUE))
#}
#st_write(dat, "dat/shp/loc.shp")

# histogram of raw data
pdf(file =  "img/hist_raw.pdf", width = 12, height = 10)
dat_clean %>% 
  dplyr::select(names(pred_lookup)) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  gather(variable, value) %>%
  left_join(., data.frame(pred_lookup, raw = names(pred_lookup)), 
    by = c("variable" = "raw")) %>%
  dplyr::select(variable = pred_lookup, value) %>% 
  mutate(variable = factor(variable, levels = pred_lookup)) %>%
  ggplot(aes(x = value)) + 
  geom_histogram(colour = "black", fill = pal[5]) + 
  facet_wrap(~variable, scales = "free") + 
  labs(x = "", y = "") + 
  theme_panel()
dev.off()

# Create bivariate relationships plot
bivar_list <- c(
  "cum_vi ~ richness",
  "cum_vi ~ evenness",
  "cum_vi ~ n_stems_gt5_ha",
  "cum_vi ~ map",
  "cum_vi ~ diurnal_temp_range",

  "s1_length ~ richness",
  "s1_length ~ evenness",
  "s1_length ~ n_stems_gt5_ha",
  "s1_length ~ map",
  "s1_length ~ diurnal_temp_range",

  "s1_green_rate ~ richness",
  "s1_green_rate ~ evenness",
  "s1_green_rate ~ n_stems_gt5_ha",
  "s1_green_rate ~ map",
  "s1_green_rate ~ diurnal_temp_range",

  "s1_senes_rate ~ richness",
  "s1_senes_rate ~ evenness",
  "s1_senes_rate ~ n_stems_gt5_ha",
  "s1_senes_rate ~ map",
  "s1_senes_rate ~ diurnal_temp_range",

  "s1_start ~ richness",
  "s1_start ~ evenness",
  "s1_start ~ n_stems_gt5_ha",
  "s1_start ~ map",
  "s1_start ~ diurnal_temp_range"
)

bivar_df <- as.data.frame(do.call(rbind, lapply(bivar_list, function(x) {
  x_var <- sym(unlist(strsplit(x, split = " ~ "))[2])
  y_var <- sym(unlist(strsplit(x, split = " ~ "))[1])

  dat_clean %>% 
    dplyr::select(!!x_var, !!y_var) %>%
    st_drop_geometry() %>%
    rename(pred = !!x_var, resp = !!y_var) %>%
    mutate(x = as.character(x_var), 
      y = as.character(y_var))
})))

bivar_df$x <- factor(bivar_df$x, levels = names(pred_lookup))
bivar_df$y <- factor(bivar_df$y, levels = names(resp_lookup))

pdf(file = "img/bivar.pdf", width = 15, height = 10)
ggplot() + 
  geom_point(data = bivar_df, aes(x = pred, y = resp), 
	colour = "black", shape = 21) +
  geom_line(data = bivar_df, aes(x = pred, y = resp),
	stat = "smooth", method = "lm", se = FALSE, colour = pal[1]) + 
  geom_line(data = bivar_df, aes(x = pred, y = resp), 
	stat = "smooth", method = "loess", colour = pal[2], se = FALSE) + 
  facet_grid(y~x, scales = "free", 
    labeller = labeller(y = resp_lookup, x = pred_lookup)) +  
  theme_panel() + 
  labs(x = "", y = "")
dev.off()

# Standardise variables
dat_std <- dat_clean %>% 
  mutate_at(.vars = c(
      "richness",
      "evenness", 
      "n_stems_gt5_ha",
      "pcoa_1", "pcoa_2", "pcoa_3", 
      "map",
      "diurnal_temp_range"),
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  dplyr::select(ends_with("_std"), names(resp_lookup), geometry) %>%
  rename_at(.vars = vars(ends_with("_std")), 
    .funs = list(~gsub("_std", "", .))) %>%
  st_transform(., UTMProj4("35S")) %>%
  mutate(x = c(unname(st_coordinates(.)[,1])),
    y = c(unname(st_coordinates(.)[,2]))) %>%
  st_drop_geometry()

# Check for collinearity
pdf(file = "img/corrplot.pdf", height = 8, width = 8)
corrPlot(dat_std[,names(pred_lookup)]) + 
  scale_x_discrete(labels = pred_lookup) + 
  scale_y_discrete(labels = pred_lookup)
dev.off()
##' None are correlated over r = 0.7, so no serious collinearity

# Linear models ----

# Define model function
phen_mod <- function(var, pre) {
  # Raw data by group
  pdf(file = paste0("img/", pre, "_richness.pdf"), height = 8, width = 10)
  bivar <- ggplot(data = dat_std, 
    aes_string(x = "richness", y = var)) + 
    geom_point() + 
    stat_smooth(method = "lm", se = TRUE)
  print(bivar)
  dev.off()

  # Define maximal model
  max_mod <- lm(get(var) ~ richness + evenness + n_stems_gt5_ha + 
    pcoa_1 + pcoa_2 + pcoa_3 + 
    map + diurnal_temp_range, 
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
    richness + evenness + n_stems_gt5_ha + 
    pcoa_1 + pcoa_2 + pcoa_3 + 
    map + diurnal_temp_range,
  data = dat_std)

  div_mod_ml <- lm(get(var) ~ richness + evenness, data = dat_std)

  env_mod_ml <- lm(get(var) ~ map + diurnal_temp_range, data = dat_std)

  mod_ml_list <- mget(c("max_mod_ml", "div_mod_ml", "env_mod_ml"))

  mod_compare_anova <- eval(parse(text=paste("anova(",
        paste("mod_ml_list[[",1:length(mod_ml_list),"]]",sep="",collapse=","),")")))

  capture.output(mod_compare_anova, file = paste0("out/", pre, "_mod_compare.txt"))
}

# Run function for key responses
phen_mod("cum_vi", "c")
phen_mod("s1_length", "l")
phen_mod("s1_green_rate", "r")
phen_mod("s1_senes_rate", "d")
phen_mod("s1_start", "s")

# spaMM models with spatial autocorrelation ----

# Fit models
max_mod_spamm_c <- fitme(cum_vi ~ richness + evenness + n_stems_gt5_ha + 
    pcoa_1 + pcoa_2 + pcoa_3 + map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_c <- fitme(cum_vi ~ map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

max_mod_spamm_l <- fitme(s1_length ~ richness + evenness + n_stems_gt5_ha + 
    pcoa_1 + pcoa_2 + pcoa_3 + map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_l <- fitme(s1_length ~ map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

max_mod_spamm_r <- fitme(s1_green_rate ~ richness + evenness + n_stems_gt5_ha + 
    pcoa_1 + pcoa_2 + pcoa_3 + map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_r <- fitme(s1_green_rate ~ map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

max_mod_spamm_d <- fitme(s1_senes_rate ~ richness + evenness + n_stems_gt5_ha + 
    pcoa_1 + pcoa_2 + pcoa_3 + map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_d <- fitme(s1_senes_rate ~ map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

max_mod_spamm_s <- fitme(s1_start ~ richness + evenness + n_stems_gt5_ha + 
    pcoa_1 + pcoa_2 + pcoa_3 + map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
null_mod_spamm_s <- fitme(s1_start ~ map + diurnal_temp_range + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")

# Define model summary function
phen_spamm_mod <- function(max_mod, null_mod, pre) {

  # Model summary
  mod_summ <- summary(max_mod)
  mod_lrt <- LRT(max_mod, null_mod)
  max_aic <- AIC(max_mod)
  max_ef_dof <- extractAIC(max_mod)
  null_aic <- AIC(null_mod)
  mod_eff <- spammEff(max_mod)

  # Hierarchical partitioning
  resp <- as.character(max_mod$predictor[[2]])

  mod_hier <- hier.part(dat_std[[resp]], dat_std %>% 
    dplyr::select(richness, evenness, n_stems_gt5_ha, pcoa_1, 
      pcoa_2, pcoa_3, map, diurnal_temp_range), barplot = FALSE)

  # Estimate degree of spatial autocorrelation
  dd <- dist(dat_std[,c("x","y")])
  mm <- MaternCorr(dd, 
    nu = max_mod$corrPars$`1`$nu, rho = max_mod$corrPars$`1`$rho)
  dist_df <- data.frame(dd = as.numeric(dd), mm = as.numeric(mm))
  dist_df <- dist_df[order(dist_df$dd),]
  pdf(file = paste0("img/", pre, "_max_mod_spamm_vario.pdf"), width = 4, height = 4)
  plot(dist_df, type = "l",
    xlab = "Pairwise distance (m)", 
    ylab = "Estimated correlation")
  dev.off()

  # Simulate residuals for model diagnostics
  sims <- simulateResiduals(max_mod)
  pdf(file = paste0("img/", pre, "_max_mod_spamm_resids.pdf"), width = 8, height = 6)
  plot(sims)
  dev.off()

  # Return list of model statistics
  ##' 1. Model predictions 
  ##' 2. Model summary 
  ##' 3. Effect sizes 
  ##' 4. Maximal model effective DoF
  ##' 5. Maximal model AIC 
  ##' 6. Null model AIC 
  ##' 7. Likelihood Ratio Test 
  ##' 8. Hierarchical partitioning
  return(list(mod_pred, mod_summ, mod_eff, max_ef_dof, 
      max_aic, null_aic, mod_lrt, mod_hier))
}

# Summarise models
spamm_c <- phen_spamm_mod(max_mod_spamm_c, null_mod_spamm_c, "c")
spamm_l <- phen_spamm_mod(max_mod_spamm_l, null_mod_spamm_l, "l")
spamm_r <- phen_spamm_mod(max_mod_spamm_r, null_mod_spamm_r, "r")
spamm_d <- phen_spamm_mod(max_mod_spamm_d, null_mod_spamm_d, "d")
spamm_s <- phen_spamm_mod(max_mod_spamm_s, null_mod_spamm_s, "s")

spamm_list <- list(spamm_c, spamm_l, spamm_r, spamm_d, spamm_s)

saveRDS(spamm_list, "dat/spamm_list.rds")
spamm_list <- readRDS("dat/spamm_list.rds")

# Export model statistics table
spamm_stat_df <- as.data.frame(do.call(rbind, lapply(spamm_list, function(x) {
      ##' If positive, max mod better
      daic_m <- x[[6]][1] - x[[5]][1]
      daic_c <- x[[6]][2] - x[[5]][2]
      ef_dof <- x[[4]][1]
  unlist(c(ef_dof, daic_m, daic_c, x[[7]]$basicLRT[1], x[[7]]$basicLRT[3]))
})))
spamm_stat_df$resp <- resp_lookup
spamm_stat_df <- spamm_stat_df[,c(length(spamm_stat_df), seq(length(spamm_stat_df) - 1))]

spamm_stat_df$p_value <- pFormat(spamm_stat_df$p_value, print_p = FALSE)

spamm_stat_tab <- xtable(spamm_stat_df, 
  label = "spamm_stat",
  align = "rrccccc",
  display = c("s", "s", "d", "d", "d", "f", "s"),
  caption = "Model fit statistics for each phenological metric.")

fileConn <- file("out/spamm_stat.tex")
writeLines(print(spamm_stat_tab, include.rownames = FALSE, 
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Export model slopes plot
spamm_eff_df <- do.call(rbind, lapply(spamm_list, function(x) {
  resp <- gsub("\\s~.*", "", as.character(x[[7]][1]$nullfit$call)[2])
  out <- x[[3]]
  out$resp <- resp
  row.names(out) <- NULL
  out[,c(4,1,2,3)]
}))

spamm_eff_df$resp <- factor(spamm_eff_df$resp, levels = names(resp_lookup))

spamm_eff_df$var <- factor(spamm_eff_df$var, levels = rev(names(pred_lookup)),
  labels = rev(pred_lookup))

pdf(file = "img/mod_spamm_slopes.pdf", width = 14, height = 6)
ggplot() + 
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_errorbarh(data = spamm_eff_df, 
    aes(xmin = est - se, xmax = est + se, y = var), height = 0.1) + 
  geom_point(data = spamm_eff_df, aes(x = est, y = var), 
    shape = 21, fill = pal[5], colour = "black") + 
  facet_wrap(~resp, scales = "free_x", nrow = 1,
    labeller = labeller(resp = resp_lookup)) +  
  theme_panel() + 
  theme(panel.grid.major.y = element_line(colour = pal[6]),
    panel.spacing = unit(2.5, "lines")) + 
  labs(x = "Slope", y = "")
dev.off()

# Export hierarchical partitioning results
mod_hier_all <- as.data.frame(do.call(cbind, lapply(spamm_list, function(x) {
  x[[8]]$I.perc[,1]
  })))
mod_hier_all$var <- pred_lookup 
mod_hier_all_clean <- mod_hier_all[,c(length(mod_hier_all),seq(length(mod_hier_all) - 1))]
names(mod_hier_all_clean) <- c("Predictor", resp_lookup)

mod_hier_tab <- xtable(mod_hier_all_clean, 
  caption = "Proportional independent contribution of each predictor in the maximal model for each phenological metric, according to hierarchical partitioning.")
align(mod_hier_tab) <- "rrccccc"

fileConn <- file("out/hier_part.tex")
writeLines(print(mod_hier_tab, include.rownames = FALSE, 
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Plot PCOA with arrows
pcoa_gather <- dat_clean %>% 
  dplyr::select(pcoa_1, pcoa_2, pcoa_3) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  gather(key, val, -pcoa_1) %>%
  mutate(key = gsub("pcoa_", "PCoA ", .$key))

##' For each pair of axes: 1,2; 1,3;
##' Calculate vector length and sort rows for 10 largest vectors
pcoa_arrows_gather <- do.call(rbind, lapply(list(c(1,2), c(1,3)), function(x) {
  # Subset pairs of axes
  pcoa_fil <- pcoa_arrows[,x]

  # Calculate vector distances
  dists <- list()
  for (i in 1:nrow(pcoa_fil)) {
    dists[i] <- as.vector(dist(t(matrix(c(pcoa_fil[i,1:2], 0,0), 2))))
  }

  # Order arrows dataframe by vector length 
  pcoa_order <- as.data.frame(pcoa_fil[order(unlist(dists), 
      decreasing = TRUE),])

  pcoa_order$sp <- gsub("\\s", "\n", row.names(pcoa_order))
  pcoa_order$key <- gsub("Axis\\.", "PCoA ", names(pcoa_order[2]))
  row.names(pcoa_order) <- NULL
  names(pcoa_order) <- c("pcoa_1", "val", "sp", "key")
  return(pcoa_order[1:10,c(3,4,1,2)])
}))

arrow_scale = 75
pdf(file = "img/pcoa.pdf", width = 10, height = 6)
ggplot() + 
  geom_point(data = pcoa_gather, 
    aes(x = pcoa_1, y = val), shape = 21, fill = pal[5], colour = "black") + 
  geom_segment(data = pcoa_arrows_gather, 
    aes(xend = pcoa_1 / arrow_scale, yend = val / arrow_scale),
    x = 0, y = 0, colour = pal[2],
    arrow = arrow(length = unit(3, "mm"))) + 
  geom_label_repel(data = pcoa_arrows_gather,
    aes(x = pcoa_1 / arrow_scale, y = val / arrow_scale, label = sp),
    label.padding = 0.085,
    segment.colour = pal[1], 
    min.segment.length = 0,
    box.padding = 0.5,
    size = 2) + 
  facet_wrap(~key) + 
  coord_equal() + 
  theme_panel() +
  labs(x = "PCoA 1", y = "")
dev.off()

# Plot location map
s1_length_tile <- as.data.frame(
  as(phen_stack$X3, "SpatialPixelsDataFrame")  # s1_length
)

pdf(file = "img/plot_loc.pdf", height = 8, width = 10)
ggplot() +
  geom_tile(data = s1_length_tile, aes(x = x, y = y, fill = X3)) +
  geom_sf(data = zambia, colour = "black", fill = NA) +
  geom_sf(data = dat_clean, colour = "black", fill = pal[6], shape = 24, size = 3) +
  scale_fill_viridis(name = "Season length\n(days)", limits = c(100, 300)) + 
  theme_panel() + 
  labs(x = "", y = "")
dev.off()

# How does species _composition_ affect phenology?

##' PCOA 1 has negative influence on season length, PCOA 2 has positive 
spamm_list[[2]][[2]]$beta_table

##' PCOA 1, 2 and 3 have positive influence on rate of greening
spamm_list[[3]][[2]]$beta_table

##' PCOA 1 has positive influence on season start
spamm_list[[5]][[2]]$beta_table

# What to export to show spatial autocorrelation effects?
##' Rho and nu from model summary objects?
