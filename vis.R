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
library(RcppRoll)
library(lubridate)

source("./plot_func.R")
source("./tex_func.R")

# Import data 
plots <- readRDS("./dat/plots.rds")
div <- readRDS("./dat/div.rds")
stat_avg <- readRDS("./dat/stat_avg.rds")
stat_all <- readRDS("./dat/stat_all.rds")
bioclim <- readRDS("./dat/bioclim.rds")
indval <- readRDS("./dat/indval.rds")
trmm <- readRDS("./dat/trmm_ts.rds")

bioclim_zambia <- readRDS("./dat/bioclim_zambia.rds")
af <- st_read("./dat/africa_countries/africa.shp")
slen_zambia <- rast("./dat/slen_zambia.tif")

# Combine plots dataframes
dat <- plots %>% 
  inner_join(., div, by = "plot_cluster") %>% 
  inner_join(., stat_avg, by = "plot_cluster") %>%
  inner_join(., bioclim, by = "plot_cluster")

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
  dplyr::select(Dormancy2, trmm_end, cluster) %>%
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

# Create facetted boxplot of MAP and MAT
mat_box <- ggplot(dat, aes(x = cluster, y = mat)) + 
  geom_boxplot(aes(fill = cluster)) + 
  scale_fill_manual(name = "Vegetation type", values = clust_pal) + 
  theme_bw() + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = NULL, y = expression("Mean Annual Temperature"~(degree*C)))

map_box <- ggplot(dat, aes(x = cluster, y = map)) + 
  geom_boxplot(aes(fill = cluster)) + 
  scale_fill_manual(name = "Vegetation type", values = clust_pal) + 
  theme_bw() + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = NULL, y = expression("Mean Annual Precipitation (mm)"))

pdf(file = "img/box_facet_map_mat.pdf", width = 8, height = 5)
mat_box + map_box + 
  plot_layout(ncol = 2)
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
    mat_mean = mean(mat, na.rm = TRUE),
    mat_sd = sd(mat, na.rm = TRUE)) %>%
  mutate(
    richness = paste0(sprintf("%.0f", richness_median), "(", sprintf("%.0f", richness_iqr), ")"),
    map = paste0(sprintf("%.0f", map_mean), "(", sprintf("%.1f", map_sd), ")"),
    mat = paste0(sprintf("%.0f", mat_mean), "(", sprintf("%.1f", mat_sd), ")")) %>%
  dplyr::select(cluster, n_sites, richness) %>%
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
  2:3] <- ""

names(clust_summ) <- c("Vegetation type", "N sites", "Richness", 
  "Indicator species", "Indicator value")

# Export indval table
clust_summ_xtable <- xtable(clust_summ,
  label = "clust_summ",
  align = rep("c", 6),
  display = rep("s", 6),
  caption = "Dufr\\^{e}ne-Legendre indicator species analysis for the vegetation type clusters identified by the PAM algorithm, based on basal area weighted species abundances \\citep{Dufrene1997}. The three species per cluster with the highest indicator values are shown alongside the median and interquartile range of site species richness and the number of sites within each cluster.")

fileConn <- file("out/clust_summ.tex")
writeLines(print(clust_summ_xtable, include.rownames = FALSE,
    table.placement = "H",
    hline.after = c(-1,0,seq(from = 3, by = 3, length.out = length(unique(indval$cluster)))),
    sanitize.text.function = function(x) {x}),
  fileConn)
close(fileConn)

# Density plots of phenological metrics per cluster
pdf(file = "img/phen_dens_clust.pdf", width = 14, height = 7)
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
  scale_colour_manual(name = "Vegetation type", values = clust_pal) + 
  theme_panel()
dev.off()
    
# Scatter plots comparing each phenological metric
bivar_plot_list <- apply(combn(names(resp_lookup), 2), 2, function(x) {
big_sig <- summary(lm(dat[[x[2]]] ~ dat[[x[1]]]))$coefficients[2,4] < 0.05
big_line <- ifelse(big_sig == TRUE, 1, 2)
veg1_sig <- summary(lm(dat[[x[2]]][dat$cluster == "Uapaca miombo"] ~ 
    dat[[x[1]]][dat$cluster == "Uapaca miombo"]))$coefficients[2,4] < 0.05
veg2_sig <- summary(lm(dat[[x[2]]][dat$cluster == "Combretaceae woodland"] ~ 
    dat[[x[1]]][dat$cluster == "Combretaceae woodland"]))$coefficients[2,4] < 0.05
veg3_sig <- summary(lm(dat[[x[2]]][dat$cluster == "Julbernardia miombo"] ~ 
    dat[[x[1]]][dat$cluster == "Julbernardia miombo"]))$coefficients[2,4] < 0.05
veg4_sig <- summary(lm(dat[[x[2]]][dat$cluster == "Cryptosepalum miombo"] ~ 
    dat[[x[1]]][dat$cluster == "Cryptosepalum miombo"]))$coefficients[2,4] < 0.05

dat$veg_line <- 1
dat$veg_line <- ifelse(veg1_sig == FALSE & dat$cluster == "Uapaca miombo", 
  2, dat$veg_line)
dat$veg_line <- ifelse(veg2_sig == FALSE & dat$cluster == "Combretaceae woodland", 
  2, dat$veg_line)
dat$veg_line <- ifelse(veg3_sig == FALSE & dat$cluster == "Julbernardia miombo", 
  2, dat$veg_line)
dat$veg_line <- ifelse(veg4_sig == FALSE & dat$cluster == "Cryptosepalum miombo", 
  2, dat$veg_line)
dat$veg_line <- as.character(dat$veg_line)

  ggplot(data = st_drop_geometry(dat), 
      aes(x = .data[[x[1]]], y = .data[[x[2]]])) + 
    geom_point(aes(fill = cluster), 
      colour = "black", shape = 21, size = 1) + 
    geom_line(aes(colour = cluster, linetype = veg_line),
      stat = "smooth", method = "lm", se = FALSE, linewidth = 1.1) + 
    geom_line(stat = "smooth", method = "lm", se = FALSE, linewidth = 1.1,
      linetype = big_line) + 
    scale_fill_manual(name = "Vegetation type", values = clust_pal) + 
    scale_colour_manual(name = "Vegetation type", 
      values = brightness(clust_pal, 0.6)) +
    scale_linetype_discrete(guide = "none") + 
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

# Extract Zambia from Africa countries
zambia <- af[af$sov_a3 == "ZMB",]

# Plot location map
pdf(file = "img/site_loc.pdf", height = 8, width = 10)
(site_loc <- ggplot() +
  geom_spatraster(data = slen_zambia) + 
  scale_fill_gradient(name = "Season length (d)", 
    low = "#1c1c1c" , high = pal[6], limits = c(100, 330), na.value = NA) + 
  new_scale_fill() +
  geom_sf(data = zambia, colour = "black", fill = NA) +
  geom_sf(data = dat, aes(fill = cluster),
    colour = "black", shape = 24, size = 2) +
  theme_panel() + 
  scale_fill_manual(name = "Vegetation type", values = clust_pal) + 
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
  scale_fill_manual(name = "Vegetation type", values = clust_pal) + 
  stat_ellipse(data = dat,
    aes(x = mat, y = map, colour = cluster), 
    type = "t", level = 0.95, linewidth = 1.2, show.legend = FALSE) + 
  scale_colour_manual(name = "Vegetation type", values = clust_pal) + 
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

# Visualise distributions of Greenup, Senescence
date_vars <- c("Greenup", "MidGreenup", "Maturity", "Senescence", 
  "MidGreendown", "Dormancy2")

pdf(file = "img/modis_distrib.pdf", width = 8, height = 12)
stat_all %>% 
  dplyr::select(plot_cluster, all_of(date_vars)) %>%
  pivot_longer(-plot_cluster) %>% 
  mutate(
    name = factor(name, levels = date_vars),
    value = as.Date(value, origin = "1970-01-01"),
    plot_num = as.numeric(as.factor(plot_cluster))) %>% 
  ggplot(., aes(x = value, y = plot_num)) + 
    geom_point(aes(colour = name)) + 
    scale_colour_scico_d(name = "Metric", palette = "roma") + 
    scale_x_date(date_labels = "%b") +
    theme_bw() +
    labs(x = "Date", y = NULL) + 
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank())
dev.off()

# Histograms of phenology metrics
pdf(file = "img/modis_hist_facet.pdf", width = 8, height = 12)
stat_all %>% 
  dplyr::select(plot_cluster, all_of(date_vars)) %>%
  pivot_longer(-plot_cluster) %>% 
  mutate(
    name = factor(name, levels = date_vars),
    value = as.Date(value, origin = "1970-01-01")) %>% 
  ggplot(., aes(x = value)) + 
    geom_histogram(aes(fill = name), colour = "black") + 
    facet_wrap(~name, scales = "fixed", ncol = 1) + 
    scale_x_date(date_labels = "%b") +
    theme_bw() +
    theme(legend.position = "none") + 
    labs(x = "Days after 1st January", y = "Frequency") 
dev.off()

# Histograms of derived metrics
stat_all %>% 
  dplyr::select(plot_cluster, names(resp_lookup)) %>% 
  pivot_longer(-plot_cluster) %>% 
  mutate(variable = factor(name, 
      levels = names(resp_lookup),
      labels = resp_lookup)) %>%
  ggplot(., aes(x = value)) + 
    geom_histogram(aes(fill = name), colour = "black") + 
    facet_wrap(~name, scales = "free") + 
    theme_bw() +
    theme(legend.position = "none") 

# Breakdown of positive senescence lag and pre-rain green-up by year
lag_table <- table(list(
    "pos green-up" = stat_all$start_lag > 0, 
    "pos senescence lag" = stat_all$end_lag > 0)) 
addmargins(lag_table)
n_table <- nrow(stat_all)

lag_table_per <- lag_table / n_table * 100

pos_sen <- round(lag_table[1,2])
pos_sen_per <- round(lag_table_per[1,2], 1)

pos_gre <- round(lag_table[2,2])
pos_gre_per <- round(lag_table_per[2,2], 1)

neg_sen <- round(sum(lag_table[,1]))
neg_sen_per <- round(sum(lag_table_per[,1]), 1)

# Compare methods of estimating start and end of rainy season
trmm_start_comp <- ggplot(stat_all, aes(x = trmm_start, y = trmm_start10)) + 
  geom_point(alpha = 0.2) +
  geom_abline(colour = "red", linetype = 2) +
  theme_bw() + 
  ggtitle("Rainy season start") + 
  labs(
    x = ">10+20 days >20+20 mm rain",
    y = "10th percentile")
# My method produces much early starts of season

trmm_end_comp <- ggplot(stat_all, aes(x = trmm_end, y = trmm_end95)) + 
  geom_point(alpha = 0.2) +
  geom_abline(colour = "red", linetype = 2) + 
  theme_bw() + 
  ggtitle("Rainy season end") + 
  labs(
    x = "<4 days 0.5 mm rain over 30 days",
    y = "95th percentile")
# My method always produces an end of season later than the q95 method

pdf(file = "./img/trmm_start_end_comp.pdf", width = 10, height = 5)
wrap_plots(trmm_start_comp, trmm_end_comp)
dev.off()

write(
  c(
    commandOutput(n_table, "nTable"),
    commandOutput(pos_sen, "posSen"),
    commandOutput(pos_sen_per, "posSenPer"),
    commandOutput(pos_gre, "posGre"),
    commandOutput(pos_gre_per, "posGrePer"),
    commandOutput(neg_sen, "negSen"),
    commandOutput(neg_sen_per, "negSenPer")
    ),
  file = "out/vis_vars.tex")
