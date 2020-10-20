# Testing effect of species composition and richness on greening with Zambia ILUAii data 
# John Godlee (johngodlee@gmail.com)
# 2020-07-29

# Packages
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(shades)
library(ggnewscale)
library(gridExtra)
library(xtable)
library(raster)
library(vegan)
library(lme4)
library(purrr)
library(MuMIn)
library(ggeffects)

source("functions.R")

# Import data 
dat <- readRDS("dat/plots_div.rds")

dat <- st_as_sf(dat)

phen_stack <- readRDS("dat/vipphen_stack.rds")

af <- st_read("/Users/johngodlee/Desktop/africa_countries/africa.shp")
zambia <- af %>% 
  filter(sov_a3 == "ZMB")

# Define variable name translation lookup
pred_lookup <- c("Richness", "Evenness", "Tree density", "MAP",
  "Wet season precip", "Dry season precip",
    "Vegetation type", "Diurnal dT")
names(pred_lookup) <- c("richness", "evenness", "n_stems_gt10_ha", "map",
  "cum_precip_seas", "cum_precip_pre",
  "cluster", "diurnal_temp_range")

resp_lookup <- c("Cumulative EVI", "Season length", 
  "Greening rate", "Senescence rate",
  "Start lag", "End lag")
names(resp_lookup) <- c("cum_vi", "s1_length", 
  "s1_green_rate", "s1_senes_rate", 
  "start_lag", "end_lag")

clust_lookup <- c("1", "2", "3", "4")
names(clust_lookup) <- c("1", "2", "3", "4")

# Calculate some statistics
dat_clean <- dat %>%
  mutate(cluster = factor(cluster, labels = clust_lookup),
    start_lag = s1_start - trmm_start,
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
  file="out/vars.tex", append=TRUE)

# Density plots of phenolgical metrics per cluster
pdf(file =  "img/phen_dens_clust.pdf", width = 12, height = 10)
dat_clean %>%
  dplyr::select(names(resp_lookup), cluster) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  gather(variable, value, -cluster) %>%
  left_join(., data.frame(resp_lookup, raw = names(resp_lookup)), 
    by = c("variable" = "raw")) %>%
  dplyr::select(variable = resp_lookup, value, cluster) %>% 
  mutate(variable = factor(variable, levels = resp_lookup)) %>%
  ggplot(., aes(x = value, colour = cluster)) + 
  geom_density(size = 1.5) + 
  facet_wrap(~variable, scales = "free") + 
  labs(x = "", y = "") +
  scale_colour_manual(values = clust_pal) + 
  theme_panel()
dev.off()

# Create bivariate relationships plot
bivar_list <- c(
  "cum_vi ~ richness",
  "cum_vi ~ evenness",
  "cum_vi ~ n_stems_gt10_ha",
  "cum_vi ~ map",
  "cum_vi ~ diurnal_temp_range",

  "s1_length ~ richness",
  "s1_length ~ evenness",
  "s1_length ~ n_stems_gt10_ha",
  "s1_length ~ map",
  "s1_length ~ diurnal_temp_range",

  "s1_green_rate ~ richness",
  "s1_green_rate ~ evenness",
  "s1_green_rate ~ n_stems_gt10_ha",
  "s1_green_rate ~ map",
  "s1_green_rate ~ diurnal_temp_range",

  "s1_senes_rate ~ richness",
  "s1_senes_rate ~ evenness",
  "s1_senes_rate ~ n_stems_gt10_ha",
  "s1_senes_rate ~ map",
  "s1_senes_rate ~ diurnal_temp_range",

  "start_lag ~ richness",
  "start_lag ~ evenness",
  "start_lag ~ n_stems_gt10_ha",
  "start_lag ~ map",
  "start_lag ~ diurnal_temp_range", 

  "end_lag ~ richness",
  "end_lag ~ evenness",
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
  st_drop_geometry()

# Check for collinearity
pdf(file = "img/corrplot.pdf", height = 8, width = 8)
corrPlot(dat_std[,names(pred_lookup)[which(names(pred_lookup) != "cluster")]]) + 
  scale_x_discrete(labels = pred_lookup) + 
  scale_y_discrete(labels = pred_lookup)
dev.off()
##' None are correlated over r = 0.7, so no serious collinearity


# Mixed models with cluster random effects ----

# Fit models
max_mod_c <- lmer(cum_vi ~ richness + evenness + cum_precip_seas + 
  diurnal_temp_range + (richness|cluster), 
  data = dat_std)
null_mod_c <- lmer(cum_vi ~ cum_precip_seas + diurnal_temp_range + 
  (1|cluster), 
  data = dat_std)

max_mod_l <- lmer(s1_length ~ richness + evenness + 
  cum_precip_seas + diurnal_temp_range + (richness|cluster), 
  data = dat_std)
null_mod_l <- lmer(s1_length ~ cum_precip_seas + diurnal_temp_range + 
  (1|cluster), 
  data = dat_std)

max_mod_r <- lmer(s1_green_rate ~ richness + evenness + cum_precip_pre + 
  (richness|cluster),
  data = dat_std)
null_mod_r <- lmer(s1_green_rate ~ cum_precip_pre + (1|cluster), 
  data = dat_std)

max_mod_d <- lmer(s1_senes_rate ~ richness + evenness + cum_precip_seas + 
  (richness|cluster), 
  data = dat_std)
null_mod_d <- lmer(s1_senes_rate ~ cum_precip_seas + 
  (1|cluster),
  data = dat_std)

max_mod_s <- lmer(start_lag ~ richness + evenness + 
  (richness|cluster), data = dat_std)
null_mod_s <- lmer(start_lag ~ (1|cluster), data = dat_std)

max_mod_e <- lmer(end_lag ~ richness + evenness +
  (richness|cluster), data = dat_std, na.action = na.omit)
null_mod_e <- lmer(end_lag ~ (1|cluster), data = dat_std)

# Define model summary function
modSumm <- function(max_mod, null_mod, pre) {

  # Predicted values
  mod_pred <- predict(max_mod)
  mod_pred_df <- data.frame(cluster = names(mod_pred), pred = as.numeric(unname(mod_pred)))

  # Model summary
  mod_summ <- summary(max_mod)
  rsq <- do.call(rbind, lapply(list(max_mod, null_mod), r.squaredGLMM))
  mod_stat_df <- data.frame(mod = c("max_mod", "null_mod"),
    aic = c(AIC(max_mod, null_mod)[,2]), 
    bic = c(BIC(max_mod, null_mod)[,2]), 
    r2m = rsq[,1], 
    r2c = rsq[,2],
    logl = c(logLik(max_mod)[1], logLik(null_mod)[1]))

  # Random effects
  rand_ef <- ranef(max_mod)[[1]]

  vars.m <- attr(rand_ef, "postVar")
  K <- dim(vars.m)[1]
  J <- dim(vars.m)[3]
  names.full <- dimnames(rand_ef)
  rand_se <- array(NA, c(J, K))
  for (j in 1:J) {
    rand_se[j, ] <- sqrt(diag(as.matrix(vars.m[, , j])))
    }
  dimnames(rand_se) <- list(names.full[[1]], names.full[[2]])

  ci_lvl <- 0.95
  ci <- 1 - ((1 - ci_lvl) / 2)

  rand_ef <- rownames_to_column(rand_ef)
  rand_se <- rownames_to_column(as.data.frame(rand_se))

  grp.names <- colnames(rand_ef)
  alabels <- rand_ef[["rowname"]]

  rand_df <- map_df(2:ncol(rand_ef), function(i) {
    out <- data.frame(estimate = rand_ef[[i]])

    # Calculate confidence intervals
    out$conf.low <- rand_ef[[i]] - (stats::qnorm(ci) * rand_se[[i]])
    out$conf.high <- rand_ef[[i]] + (stats::qnorm(ci) * rand_se[[i]])
    out$se <- (stats::qnorm(ci) * rand_se[[i]])

    # set column names (variable / coefficient name)
    # as group indicator, and save axis labels and title in variable
    out$facet <- grp.names[i]
    out$term <- factor(alabels)
    out$resp <- names(max_mod@frame)[1] 

    # create default grouping, depending on the effect:
    # split positive and negative associations with outcome
    # into different groups
    out$group <- dplyr::if_else(out$estimate > 0, "pos", "neg")

    return(out)
  })

  # Return list of model statistics
  ##' 1. Max mod object
  ##' 2. Null mod object
  ##' 3. Model predicted values
  ##' 4. Model summary 
  ##' 5. Model fit statistics
  ##' 6. Effect sizes for each random effect
  return(list(max_mod, null_mod, mod_pred, mod_summ, mod_stat_df, rand_df))
}

# Summarise models
summ_c <- modSumm(max_mod_c, null_mod_c, "c")
summ_l <- modSumm(max_mod_l, null_mod_l, "l")
summ_r <- modSumm(max_mod_r, null_mod_r, "r")
summ_d <- modSumm(max_mod_d, null_mod_d, "d")
summ_s <- modSumm(max_mod_s, null_mod_s, "s")
summ_e <- modSumm(max_mod_e, null_mod_e, "e")

summ_list <- list(summ_c, summ_l, summ_r, summ_d, summ_s, summ_e)

saveRDS(summ_list, "dat/summ_list.rds")
summ_list <- readRDS("dat/summ_list.rds")

# Export model statistics table
mod_stat_df <- as.data.frame(do.call(rbind, lapply(summ_list, function(x) {
      ##' If positive, max mod better
      daic <- x[[5]][2,2] - x[[5]][1,2]
      dbic <- x[[5]][2,3] - x[[5]][1,3]
      r2m <- x[[5]][1,4]
      r2c <- x[[5]][1,5]
      logl <- x[[5]][1,6]

  unlist(c(daic, dbic, r2m, r2c, logl))
}))) 
names(mod_stat_df) <- c("daic", "dbic", "r2m", "r2c", "logl")
mod_stat_df$resp <- unlist(lapply(summ_list, function(x) {  names(x[[1]]@frame)[1] }))
mod_stat_df <- mod_stat_df[,c(length(mod_stat_df), seq(length(mod_stat_df) - 1))]
mod_stat_df$resp <- gsub("_", ".", mod_stat_df$resp)
mod_stat_tab <- xtable(mod_stat_df, 
  label = "mod_stat",
  align = "rrccccc",
  display = c("s", "s", "f", "f", "f", "f", "f"),
  caption = "Model fit statistics for each phenological metric.")

fileConn <- file("out/mod_stat.tex")
writeLines(print(mod_stat_tab, include.rownames = FALSE, 
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Export model slopes plot
mod_eff_df <- do.call(rbind, lapply(summ_list, function(x) {
  resp <- names(x[[1]]@frame)[1] 
  out <- as.data.frame(x[[4]]$coefficients)
  out$resp <- resp
  out <- rownames_to_column(out, "fixeff")
  out[,c(5,1,2,3,4)]
})) %>%
  mutate(resp = factor(resp, levels = names(resp_lookup)),
    fixeff = factor(fixeff, levels = rev(names(pred_lookup)), labels = rev(pred_lookup)),
    pos_neg = if_else(Estimate >= 0, "pos", "neg")) %>%
  filter(!is.na(fixeff)) %>%
  rename(est = Estimate, se = `Std. Error`, tval = `t value`)

pdf(file = "img/mod_slopes.pdf", width = 12, height = 6)
ggplot() + 
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_errorbarh(data = mod_eff_df, 
    aes(xmin = est - se, xmax = est + se, y = fixeff, colour = pos_neg), 
    height = 0.1) + 
  geom_point(data = mod_eff_df, 
    aes(x = est, y = fixeff, fill = pos_neg), 
    shape = 21, colour = "black") + 
  scale_colour_manual(values = pal[3:4]) + 
  scale_fill_manual(values = pal[3:4]) + 
  facet_wrap(~resp, scales = "free_x", nrow = 2,
    labeller = labeller(resp = resp_lookup)) +  
  theme_panel() + 
  theme(panel.grid.major.y = element_line(colour = pal[6]),
    panel.spacing = unit(2.5, "lines"),
    legend.position = "none") + 
  labs(x = "Slope", y = "")
dev.off()

# Export random effects plots
rand_marg_df <- do.call(rbind, lapply(summ_list, function(i) {
  preds <- as.data.frame(ggpredict(model = i[[1]], 
      terms = c("richness", "cluster"), type = "random"))
  preds$resp <- names(i[[1]]@frame)[1]
  names(preds) <- c("richness", "pred", "cluster", "resp")
  return(preds)
}))

pdf(file = "img/mod_marg.pdf", width = 12, height = 6)
ggplot() + 
  geom_line(data = rand_marg_df, 
    aes(x = richness, y = pred, colour = cluster),
    size = 1.5) + 
  facet_wrap(~resp, scales = "free_y", 
    labeller = labeller(resp = resp_lookup)) +  
  scale_colour_manual(values = clust_pal) +
  theme_panel() +
  labs(x = "Species richness", y = "")
dev.off()

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
  scale_fill_manual(name = "", values = clust_pal) + 
  labs(x = "", y = "")
dev.off()


