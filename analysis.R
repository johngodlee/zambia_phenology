# Pre-processing analysis
# John Godlee (johngodlee@gmail.com)
# 2020-07-29

# Packages
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggfortify)
library(ggnewscale)
library(shades)
library(gridExtra)
library(raster)
library(xtable)

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
  mutate(
    cluster = factor(cluster, 
      labels = clust_lookup[1:length(unique(dat$cluster))]),
    start_lag = -(s1_start - trmm_start),
    end_lag = s1_end - trmm_end,
    cum_vi = cum_vi / 100000)

dat_nogeom <- st_drop_geometry(dat_clean) %>%
  filter(plot_cluster != "ZIS_1232")
  
# Create density distributions of season start and end dates
start_dens_plot <- dat_nogeom %>%
  dplyr::select(s1_start, trmm_start, cluster) %>%
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

end_dens_plot <- dat_nogeom %>%
  dplyr::select(s1_end, trmm_end, cluster) %>%
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

# Density plots of phenological metrics per cluster
pdf(file = "img/phen_dens_clust.pdf", width = 12, height = 10)
dat_nogeom %>%
  dplyr::select(names(resp_lookup), cluster) %>%
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

# PCA of phenological metrics by plot and community
phen_pca <- prcomp(dat_nogeom[,names(resp_lookup)], 
  center = TRUE, scale. = TRUE)

phen_pca_tidy <- as.data.frame(phen_pca$x) %>%
  bind_cols(dat_nogeom)

phen_pca_hulls <- findHull(phen_pca_tidy, "PC1", "PC2", "cluster")

phen_pca_arrows <- data.frame(x = rownames(phen_pca$rotation), 
  phen_pca$rotation)

pdf(file = "img/phen_pca_clust.pdf", width = 8, height = 6)
ggplot() + 
  geom_polygon(data = phen_pca_hulls, aes(x = PC1, y = PC2, colour = cluster), 
    fill = NA) + 
  geom_point(data = phen_pca_tidy, aes(x = PC1, y = PC2, fill = cluster),
    shape = 21, colour = "black") + 
  geom_segment(data = phen_pca_arrows, 
    aes(x = 0, y = 0, xend = (PC1 * 5), yend = (PC2 * 5)), 
      arrow = arrow(length = unit(1/2, "picas")), colour = "black") + 
  geom_label(data = phen_pca_arrows, 
    aes(x = (PC1 * 4), y = (PC2 * 4), label = x)) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  theme_panel()
dev.off()

# Create bivariate relationships plot
bivar_list <- c(
  "cum_vi ~ eff_rich",
  "cum_vi ~ n_stems_ge10_ha",
  "cum_vi ~ map",
  "cum_vi ~ diurnal_temp_range",
  "cum_vi ~ diam_quad_mean",

  "s1_length ~ eff_rich",
  "s1_length ~ n_stems_ge10_ha",
  "s1_length ~ map",
  "s1_length ~ diurnal_temp_range",
  "s1_length ~ diam_quad_mean",

  "s1_green_rate ~ eff_rich",
  "s1_green_rate ~ n_stems_ge10_ha",
  "s1_green_rate ~ map",
  "s1_green_rate ~ diurnal_temp_range",
  "s1_green_rate ~ diam_quad_mean",

  "s1_senes_rate ~ eff_rich",
  "s1_senes_rate ~ n_stems_ge10_ha",
  "s1_senes_rate ~ map",
  "s1_senes_rate ~ diurnal_temp_range",
  "s1_senes_rate ~ diam_quad_mean",

  "start_lag ~ eff_rich",
  "start_lag ~ n_stems_ge10_ha",
  "start_lag ~ map",
  "start_lag ~ diurnal_temp_range", 
  "start_lag ~ diam_quad_mean",

  "end_lag ~ eff_rich",
  "end_lag ~ n_stems_ge10_ha",
  "end_lag ~ map",
  "end_lag ~ diurnal_temp_range",
  "end_lag ~ diam_quad_mean"
)

bivar_df <- as.data.frame(do.call(rbind, lapply(bivar_list, function(x) {
  x_var <- sym(unlist(strsplit(x, split = " ~ "))[2])
  y_var <- sym(unlist(strsplit(x, split = " ~ "))[1])

  dat_nogeom %>% 
    dplyr::select(!!x_var, !!y_var, cluster) %>%
    rename(pred = !!x_var, resp = !!y_var) %>%
    mutate(x = as.character(x_var), 
      y = as.character(y_var))
})))

bivar_df$x <- factor(bivar_df$x, levels = names(pred_lookup))
bivar_df$y <- factor(bivar_df$y, levels = names(resp_lookup))

#bivar_df %>% 
#  filter(x == "diam_cov", y == "cum_vi") %>%
#  ggplot(., aes(x = pred, y = resp)) + 
#  geom_point()

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

# Scatter plots comparing each phenological metric
phen_bivar_df <- do.call(rbind, lapply(combn(names(resp_lookup), 2, simplify = FALSE),
  function(x) { 
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

pdf(file = "img/phen_bivar.pdf", width = 15, height = 10)
ggplot() + 
  geom_point(data = phen_bivar_df, aes(x = pred, y = resp, fill = cluster),
    colour = "black", shape = 21) + 
  geom_line(data = phen_bivar_df, aes(x = pred, y = resp),
	stat = "smooth", method = "lm", colour = "black", se = FALSE, size = 1.5) + 
  geom_line(data = phen_bivar_df, aes(x = pred, y = resp, colour = cluster), 
	stat = "smooth", method = "lm", se = FALSE, size = 1.5) + 
  facet_wrap(~pair_label, scales = "free") + 
  scale_fill_manual(name = "Cluster", values = clust_pal) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) + 
  theme_panel() + 
  labs(x = "", y = "")
dev.off()

# Correlation table comparing each phenological metric across veg. types

phen_bivar_split <- split(phen_bivar_df, phen_bivar_df$pair_name)

phen_corr_list <- unlist(lapply(phen_bivar_split, function(x) {
  clust_split <- split(x, x$cluster)

  all_corr <- cor.test(x$pred, x$resp, method = "pearson") 

  clust_corr <- lapply(clust_split, function(y) {
    cor.test(y$pred, y$resp, method = "pearson") 
  })

  out <- c(clust_corr, list(all_corr))

  names(out) <- c(as.character(seq_len(length(out)-1)), "all")

  return(out)
}), recursive = FALSE)

phen_corr_df <- do.call(rbind, lapply(seq(length(phen_corr_list)), function(x) {
  nam <- unlist(strsplit(names(phen_corr_list[x]), "-|\\."))
  corr <- phen_corr_list[[x]]$estimate
  pval <- phen_corr_list[[x]]$p.value
  dof <- phen_corr_list[[x]]$parameter
  tstat <- phen_corr_list[[x]]$statistic
  out <- data.frame(pred = nam[1], resp = nam[2], cluster = nam[3], 
    dof, corr, tstat, pval) 

  return(out)
}))

phen_corr_df$pred_pretty <- resp_lookup[match(phen_corr_df$pred, names(resp_lookup))]
phen_corr_df$resp_pretty <- resp_lookup[match(phen_corr_df$resp, names(resp_lookup))]
phen_corr_df$pval_pretty <- pFormat(phen_corr_df$pval)

phen_corr_df_out <- phen_corr_df[,c("pred_pretty", "resp_pretty", "cluster", 
  "dof", "corr", "tstat", "pval_pretty")]

names(phen_corr_df_out) <- c("X", "Y", "Cluster", "DoF", 
  "$r$", "$t$", "Prob.")

phen_corr_xtable <- xtable(phen_corr_df_out, 
  label = "phen_corr",
  align = c("c", "r", "r", "c", "r", "r", "r", "r"),
  display = c("s", "s", "s", "d", "d", "f", "f", "f"),
  caption = "Pearson's Correlation Coefficient test for all pairwise combinations of the six phenological metrics used in this study.")

nclust <- length(unique(phen_corr_df_out$Cluster))
fileConn <- file("out/phen_corr.tex")
writeLines(print(phen_corr_xtable, include.rownames = FALSE,
    table.placement = "H",
    hline.after = c(-1, 0, seq(from = nclust , by = nclust, length.out = nrow(phen_corr_df) / nclust)),
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Standardise variables
dat_std <- dat_clean %>% 
  mutate_at(.vars = names(pred_lookup)[which(names(pred_lookup) != "cluster")],
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  dplyr::select(plot_cluster, ends_with("_std"), cluster, names(resp_lookup), 
    geometry) %>%
  rename_at(.vars = vars(ends_with("_std")), 
    .funs = list(~gsub("_std", "", .))) %>%
  st_transform(., UTMProj4("35S")) %>%
  mutate(x = c(unname(st_coordinates(.)[,1])),
    y = c(unname(st_coordinates(.)[,2]))) %>%
  st_drop_geometry() %>%
  drop_na()

# Check dataframes have same rows
stopifnot(nrow(dat_std) == nrow(dat_clean))

# Check for collinearity
pdf(file = "img/corrplot.pdf", height = 8, width = 10)
corrPlot(dat_std[,names(pred_lookup)[which(names(pred_lookup) != "cluster")]]) + 
  scale_x_discrete(labels = pred_lookup) + 
  scale_y_discrete(labels = pred_lookup)
dev.off()
##' None are correlated over r = 0.7, so no serious collinearity

# MANOVA to show variation within vegetation clusters vs. outside
dat_manova <- dat_std %>% 
  mutate_at(.vars = names(resp_lookup),
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  dplyr::select(cluster, ends_with("_std"))

phen_manova <- manova(as.matrix(dat_manova[,grepl("_std", names(dat_manova))]) ~ 
  dat_manova$cluster)

# Tukey's tests for each phenological metric per cluster
resps <- names(dat_manova)[grepl("_std", names(dat_manova))]

resps_out <- do.call(rbind, lapply(resps, function(x) {
  mod <- aov(dat_manova[[x]] ~ dat_manova$cluster)
  tukey <- TukeyHSD(mod)
  tukey_clean <- rownames_to_column(as.data.frame(tukey[[1]]), "comp")
  tukey_clean$resp <- x
  return(tukey_clean)
  }))

# Are all intervals overlapping 0?
all(resps_out$upr - resps_out$lwr >= 0)

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

# How many sites are there?
write(
  commandOutput(nrow(dat_nogeom), "nSites"),
  file = "out/analysis_vars.tex")

# Save data ready for models
saveRDS(dat_std, "dat/plots_anal.rds")

