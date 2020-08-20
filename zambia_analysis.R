# Testing effect of species composition and richness on greening with Zambia ILUAii data 
# John Godlee (johngodlee@gmail.com)
# 2020-07-29

# Preamble ----

# Reset env. 
# rm(list = ls())
# dev.off()

# Delete variables TeX file 
if (file.exists("out/z_vars.tex")) {
  file.remove("out/z_vars.tex")
}

# Set working directory

# Packages
library(ggplot2)
library(ggrepel)
library(viridis)
library(dplyr)
library(tidyr)
library(sf)
library(raster)
library(gdalUtils)
library(vegan)
library(ape)
library(lme4)
library(gridExtra)
library(sjPlot)
library(effects)
library(MuMIn)
library(car)
library(hier.part)
library(DHARMa)
library(spaMM)

#' Get valid UTM zone from latitude and longitude in WGS84 decimal degrees
#'
#' @param x vector of longitude coordinates in decimal degrees
#' @param y vector of latitude coordinate in decimal degrees
#'
#' @return Vector of UTM zones for each latitude-longitude pair
#' 
#' @export
#' 
latLong2UTM <- function(x, y) {
  unlist(lapply(1:length(x), function(z) {
    paste((floor((x[z] + 180) / 6) %% 60) + 1,
      ifelse(y[z] < 0, "S", "N"),
      sep = "")
  }))
}

#' Generate a valid UTM WGS84 proj4string given a UTM zone
#'
#' @param x character vector defining UTM zones
#'
#' @return UTM proj4string character vector
#' 
#' @export
#' 
UTMProj4 <- function(x){
  unlist(lapply(1:length(x), function(y) {
    paste0(
      "+proj=utm +zone=",
      gsub("[A-z]", "", as.character(x[y])),
      ifelse(gsub("[0-9]", "", as.character(x[y])) == "S", " +south", ""),
      " +ellps=WGS84")
  }))
}

pcoa_arrows <- function(given_pcoa, trait_df, sort = FALSE) {
    n <- nrow(trait_df)
    points.stand <- scale(given_pcoa$vectors)
    
    # Compute covariance of variables with all axes
    S <- cov(trait_df, points.stand)
    
    # Select only positive eigenvalues
    pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
    
    # Standardize value of covariance (see Legendre & Legendre 1998)
    U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
    colnames(U) <- colnames(given_pcoa$vectors)
    
    # Add values of covariances inside object
    given_pcoa$U <- U

    if (sort) {
      # Compute arrow vector lengths
      dists <- list()
      for (i in 1:nrow(given_pcoa$U)) {
        dists[i] <- as.vector(dist(t(matrix(c(given_pcoa$U[i,1:2], 0,0), 2))))
      }

      # Order arrows dataframe by vector length 
      given_pcoa$U <- as.data.frame(given_pcoa$U[order(unlist(dists), 
          decreasing = TRUE),])
    }
    
    return(given_pcoa)
}

num_format <- function(x, digits = 2, method = "round"){
  sprintf(paste0("%.",digits,"f"),
    if(method == "round"){
      round(x, digits = digits)
    }else if(method == "signif"){
      signif(x, digits = digits)
    })
}

command_output <- function(x, name){ 
  paste0("\\newcommand{\\",
    ifelse(missing(name), deparse(substitute(x)), name), 
    "}{", 
    x, 
    "}"
  )
}

corrplot <- function(x, col = c("blue", "white", "red"), ...) {
  corr <- psych::corr.test(x, ...) 
  corr_ci <- data.frame(raw.lower = corr$ci$lower, raw.r = corr$ci$r, 
    raw.upper = corr$ci$upper, raw.p = corr$ci$p, 
    adj.lower = corr$ci.adj$lower.adj, adj.upper = corr$ci.adj$upper.adj)
  corr_ci$vars <- row.names(corr_ci)
  corr_ci$conf_x <- unlist(sapply(1:(length(x)-1), function(i){
      c(1:(length(x)-1))[i:(length(x)-1)]
    })) + 1
  rev_mat <- (length(x)-1):1
  corr_ci$conf_y <- unlist(sapply(1:(length(x)-1), function(i){
      rep(i, times = rev_mat[i])
    }))
  n_seq <- 2:length(x)
  corr_ci$y_var <- unlist(sapply(1:(length(x)-1), function(i){
      rep(row.names(corr[[1]])[i], rev_mat[i])
    }))
  corr_ci$x_var <- unlist(sapply(1:(length(x)-1), function(i){
      row.names(corr[[1]])[n_seq[i]:length(x)]
    }))
  corr_ci$x_var <- factor(corr_ci$x_var, levels = unique(corr_ci$x_var))
  corr_ci$y_var <- factor(corr_ci$y_var, levels = unique(corr_ci$y_var))
  corr_ci$conf <- (corr_ci$raw.lower > 0) == (corr_ci$raw.upper > 0)
  corr_ci$raw.r <- round(corr_ci$raw.r, 2)

  ggplot2::ggplot() + 
    ggplot2::geom_tile(data = corr_ci, 
      ggplot2::aes(x = x_var, y = y_var, 
        fill = raw.r), colour = "black") + 
    ggplot2::geom_text(data = corr_ci, 
      ggplot2::aes(x = x_var, y = y_var, label = raw.r),
      size = 3) + 
    ggplot2::geom_point(data = corr_ci[corr_ci$conf == FALSE,], 
      ggplot2::aes(x = x_var, y = y_var), 
      fill = NA, colour = "black", shape = 21, size = 11) + 
    ggplot2::scale_fill_gradient2(name = "r", 
      low = col[1], mid = col[2], high = col[3]) + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x = "", y = "") + 
    ggplot2::coord_equal() +
    ggplot2::theme(legend.position = "none")
}

# Import data ----

plots <- read.csv("~/git_proj/seosaw_data/data_out/plots_v2.7.csv")
stems <- read.csv("~/git_proj/seosaw_data/data_out/stems_latest_v2.7.csv")
af <- st_read("/Volumes/john/africa_countries/africa.shp")

zambia <- af %>% 
  filter(sov_a3 == "ZMB")

map_stack <- stack(list.files("/Volumes/john/worldclim/wc2.1_30s_prec", "*.tif", 
  full.names = TRUE))
#map_zambia <- crop(map_stack, zambia)
#map_sum <- calc(map_zambia, sum)
#map <- mask(map_sum, zambia)
#saveRDS(map, "data/map.rds")
map <- readRDS("data/map.rds")

# List of vipphen layers
vipphen_file_list <- list.files(path = "/Volumes/john/modis_vipphen", 
  pattern = ".*tif$", full.names = TRUE)

# Split by variable
var_list <- as.numeric(gsub("\\.tif", "", gsub(".*_", "", vipphen_file_list)))
vipphen_file_list_split <- split(vipphen_file_list, var_list)

# Stack by variable
phen_stack_list <- lapply(vipphen_file_list_split, stack)

# Select variables I care about
phen_stack_list_fil <- phen_stack_list[c(1,2,3,5,6,22,23,24,25)]

# Subset plots and data fields ----
plots_fil_sf <- plots %>% 
  filter(
    prinv == "Siampale A.") %>% # Zambia only 
  dplyr::select(
    plot_id, plot_cluster, longitude_of_centre, latitude_of_centre, 
    richness, shannon, simpson, evenness, n_stems_gt5_ha, ba_ha, agb_ha,
    mat = bio1, diurnal_temp_range = bio2, map = bio12, 
    clay = CLYPPT, sand = SNDPPT, cec = CECSOL) %>%
  group_by(plot_cluster) %>%
  summarise(
    plot_id = paste0(plot_id, collapse = ","),
    longitude_of_centre = mean(longitude_of_centre, na.rm = TRUE),
    latitude_of_centre = mean(latitude_of_centre, na.rm = TRUE),
    n_stems_gt5_ha = mean(n_stems_gt5_ha, na.rm = TRUE),
    ba_ha = mean(ba_ha, na.rm = TRUE),
    agb_ha = mean(agb_ha, na.rm = TRUE),
    mat = mean(mat, na.rm = TRUE),
    diurnal_temp_range = mean(diurnal_temp_range, na.rm = TRUE),
    map = mean(map, na.rm = TRUE),
    clay = mean(clay, na.rm = TRUE),
    sand = mean(sand, na.rm = TRUE),
    cec = mean(cec, na.rm = TRUE)) %>%
  st_as_sf(., coords = c("longitude_of_centre", "latitude_of_centre")) %>%
  `st_crs<-`(4326) %>%
  mutate(plot_id_vec = strsplit(as.character(plot_id), 
  split = ",")) %>%
  mutate(plot_id_length = sapply(.$plot_id_vec, length)) %>%
  filter(plot_id_length == 4) %>%
  dplyr::select(-plot_id_length)

# Crease plot ID / plot Cluster lookup table
plot_id_lookup <- plots %>% 
  filter(plot_id %in% unlist(plots_fil_sf$plot_id_vec)) %>%
  dplyr::select(plot_cluster, plot_id)

# Filter stems by plots
stems_fil <- stems %>%
  inner_join(., plot_id_lookup, by = "plot_id") %>%
  dplyr::select(-plot_id)

# Mopane stem percentage filter
mopane_per <- 0.5
stems_ha <- 50

write(
  c(
    command_output(stems_ha, "stemsHa"),
    command_output(mopane_per, "mopanePer")
  ),
  file="out/z_vars.tex", append=TRUE)

# Create tree species abundance matrix
tree_ab_mat <- stems_fil %>% 
  dplyr::select(plot_cluster, tree_id, species_name_clean) %>%
  filter(!is.na(species_name_clean)) %>%
  group_by(plot_cluster, tree_id) %>%
  filter(row_number() == 1) %>%
  group_by(plot_cluster, species_name_clean, .drop = FALSE) %>%
  tally() %>%
  spread(species_name_clean, n, fill = 0) %>%
  ungroup() %>%
  mutate_at(vars(-plot_cluster), as.double) %>%
  as.data.frame() %>%
  dplyr::select(-`Indet indet`)
rownames(tree_ab_mat) <- tree_ab_mat[["plot_cluster"]]
tree_ab_mat <- tree_ab_mat[,-1, drop=FALSE]
tree_ab_mat <- filter(tree_ab_mat, (tree_ab_mat$`Colophospermum mopane` / rowSums(tree_ab_mat)) < mopane_per)
tree_ab_mat <- filter_all(tree_ab_mat, any_vars(. != 0))
tree_ab_mat <- tree_ab_mat[stems_ha < rowSums(tree_ab_mat) / 0.4,]

# Remove plots not in tree abundance matrix
plots_fil_sf <- filter(plots_fil_sf, plot_cluster %in% rownames(tree_ab_mat))

# Write as shapefile
if (file.exists("data/z_shp/z_loc.shp")) {
  file.remove(list.files("data/z_shp", "z_loc.*", full.names = TRUE))
}
st_write(plots_fil_sf, "data/z_shp/z_loc.shp")

write(
  command_output(nrow(plots_fil_sf), "nSites"),
  file="out/z_vars.tex", append=TRUE)

# Create plot extent
plots_extent <- extent(plots_fil_sf) 
plots_extent[c(1,3)] <- plots_extent[c(1,3)] -1
plots_extent[c(2,4)] <- plots_extent[c(2,4)] +1

plots_bbox <- st_as_sfc(st_bbox(plots_extent)) %>%
  `st_crs<-`(4326)

pdf(file = "img/z_loc_map.pdf", height = 10, width = 8)
ggplot() + 
  geom_sf(data = af) + 
  geom_sf(data = plots_bbox, fill = "blue", alpha = 0.5) + 
  geom_sf(data = plots_fil_sf) 
dev.off()

# Add VIPPHEN data to plots ----

# Crop phenology
phen_crop_list <- lapply(phen_stack_list_fil, crop, y = plots_extent)

# Set CRS of raster to WGS84
phen_wgs_list <- lapply(phen_crop_list, projectRaster, 
  crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

# Save image of VIPPHEN
phen_zam_list <- lapply(phen_wgs_list, mask, zambia)
phen_zam_mean_stack <- stack(lapply(phen_zam_list, calc, mean))

s1_length_tile <- as.data.frame(
  as(phen_zam_mean_stack$X3, "SpatialPixelsDataFrame")  # s1_length
)

pdf(file = "img/z_plot_loc.pdf", height = 10, width = 13)
ggplot() +
  geom_tile(data = s1_length_tile, aes(x = x, y = y, fill = X3)) +
  geom_sf(data = zambia, colour = "black", fill = NA) +
  geom_sf(data = plots_fil_sf, colour = "black", fill = "white", shape = 24) +
  scale_fill_viridis(name = "Growth season\nlength (days)", limits = c(100, 300)) + 
  theme_classic() + 
  labs(x = "", y = "")
dev.off()

# Pull values and take mean across years
plots_phen_extract <- as.data.frame(do.call(cbind, lapply(phen_wgs_list, function(x) { raster::extract(x, plots_fil_sf, method = "simple") %>%
  as.data.frame(.) %>%
  mutate_all(~ case_when(
      . == -1 ~ NA_real_,
      . < -5000 ~ NA_real_,
      TRUE ~ .)) %>%
  rowMeans(., na.rm = TRUE)
})))

names(plots_phen_extract) <- c("s1_start", "s1_end", "s1_length", 
  "s1_green_rate", "s1_senes_rate", "cum_vi", "avg_vi", "bg_vi", "n_seasons")

plots_phen <- plots_fil_sf %>% 
  bind_cols(., plots_phen_extract) %>%
  filter(n_seasons < 2) %>%
  dplyr::select(-n_seasons) %>%
  filter(!is.na(s1_length))

# Save intermediate as .rds ----
saveRDS(plots_phen, file = "data/z_phen.rds")
plots_phen <- readRDS("data/z_phen.rds")

# Get diversity data ----
# Run PCOA on distance matrix
tree_dist <- vegdist(tree_ab_mat, method = "bray")
tree_pcoa <- pcoa(tree_dist)

# Percentage variance explained
eig_per <- list()
for (i in 1:length(tree_pcoa$values$Relative_eig)) {
  eig_per[i] <- sum(tree_pcoa$values$Relative_eig[1:i]) / sum(tree_pcoa$values$Relative_eig)
}

# First three axes
write(
  command_output(num_format(sum(unlist(eig_per[1:3]))), "pcoaPer"),
  file="out/z_vars.tex", append=TRUE)

# Extract arrow vectors for species
trait_pcoa_arrows <- pcoa_arrows(tree_pcoa, tree_ab_mat, sort = TRUE)

arrows_main <- trait_pcoa_arrows$U
arrows_main$species <- row.names(arrows_main)
arrows_main$species <- gsub(" ", "\n", arrows_main$species)

# Extract values
pcoa_values <- as.data.frame(tree_pcoa$vectors)[,1:5]
names(pcoa_values) <- paste("pcoa", 1:5, sep = "_")
pcoa_values$plot_cluster <- row.names(pcoa_values)

# Calculate common  diversity statistics
div_df <- data.frame(plot_cluster = row.names(tree_ab_mat), 
  richness = unname(rowSums(tree_ab_mat != 0)),
  shannon = diversity(tree_ab_mat),
  simpson = diversity(tree_ab_mat, "simpson"),
  evenness = diversity(tree_ab_mat, "invsimpson") / rowSums(tree_ab_mat > 0))

# Add values to data
plots_div <- plots_phen %>%
  left_join(., pcoa_values, by = "plot_cluster") %>%
  left_join(., div_df, by = "plot_cluster") %>% 
  filter(!is.na(richness))

# Plot PCOA with arrows
pdf(file = "img/z_pcoa_arrows.pdf", width = 10, height = 8)
ggplot() +
  geom_point(data = plots_div, aes(x = pcoa_1, y = pcoa_2, colour = s1_length)) +
  geom_segment(data = arrows_main[1:10,],
    aes(xend = Axis.1 / 75, yend = Axis.2 / 75),
    x = 0, y = 0, alpha = 0.7, colour = "red",
    arrow = arrow(length = unit(3, "mm"))) + 
  geom_label_repel(data = arrows_main[1:10,],
    aes(x = Axis.1 / 75, y = Axis.2 / 75, label = species),
    label.padding = unit(0.1, "lines"), size = 3) + 
  scale_colour_viridis(name = "Growth season\nlength (days)") +
  theme_classic() + 
  labs(x = "PCo 1", y = "PCo 2")
dev.off()

# Plot PCOA axes
pcoa_gather <- plots_div %>%
  dplyr::select(pcoa_1, pcoa_2, pcoa_3, pcoa_4, pcoa_5, s1_length) %>%
  gather(key, val, -pcoa_1, -s1_length) %>%
  mutate(key = toupper(gsub("oa_", "o ", .$key)))

pdf(file = "img/z_pcoa.pdf", width = 10, height = 8)
ggplot() + 
  geom_point(data = pcoa_gather, 
    aes(x = pcoa_1, y = val, colour = s1_length)) + 
  scale_colour_viridis(name = "Growth season\nlength (days)") + 
  facet_wrap(~key) + 
  coord_equal() + 
  theme_bw() + 
  labs(x = "PCo 1", y = "")
dev.off()

# Save data as .rds
saveRDS(plots_div, "data/z_div.rds")

# Analysis ----

dat <- readRDS("data/z_div.rds")

pdf(file =  "img/z_hist_raw.pdf", width = 12, height = 10)
dat %>% 
  dplyr::select(-geometry, -plot_id, -plot_id_vec, -plot_cluster) %>%
  gather(variable, value) %>%
  ggplot(aes(x = value)) + 
  geom_histogram(colour = "black", fill = "grey") + 
  facet_wrap(~variable, scales = "free") + 
  labs(x = "", y = "") + 
  theme_bw()
dev.off()

# Explore data with bivariate plots ----

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
  lm(eval(x), data = dat)
})

# Create plots
plot_list <- lapply(bivar_list, function(x) {
  xvar <- unlist(strsplit(x, split = " ~ "))[2]
  yvar <- unlist(strsplit(x, split = " ~ "))[1]
  
  ggplot() + 
    geom_point(data = dat,
	  aes_string(x = xvar, y = yvar),
	  colour = "black", shape = 21, alpha = 0.8) + 
	geom_line(data = dat,
	  aes_string(x = xvar, y = yvar),
	  stat = "smooth", method = "lm", se = FALSE, colour = "blue") + 
	geom_line(data = dat,
	  aes_string(x = xvar, y = yvar), 
	  stat = "smooth", method = "loess", colour = "#DE6400", se = FALSE) + 
	theme_bw()
})

# Arrange on grid
n <- length(plot_list)
n_col <- floor(sqrt(n))

pdf(file =  "img/z_bivar.pdf", width = 12, height = 10)
do.call("grid.arrange", c(plot_list, ncol = 6))
dev.off()

dat_std <- dat %>% 
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
  dplyr::select(ends_with("_std"), s1_length, cum_vi, bg_vi,
    s1_green_rate, s1_senes_rate, s1_start, geometry) %>%
  st_transform(., UTMProj4("35S")) %>%
  mutate(x = c(unname(st_coordinates(.)[,1])),
    y = c(unname(st_coordinates(.)[,2]))) %>%
  st_drop_geometry()

# Check for collinearity
pdf(file = "img/z_corrplot.pdf", height = 8, width = 8)
corrplot(dat_std[,c("richness_std", "evenness_std", "n_stems_gt5_ha_std", 
  "pcoa_1_std", "pcoa_2_std", "pcoa_3_std",
  "map_std", "diurnal_temp_range_std")])
dev.off()
##' None are correlated over r = 0.7, so no serious collinearity

pdf(file = "img/corrplot_response.pdf", height = 8, width = 8)
corrplot(dat_std[,c("s1_length", "s1_green_rate", "s1_senes_rate", "cum_vi", "bg_vi")])
dev.off()

# Linear models ----

# TESTING USING spaMM models
max_mod <- fitme(s1_length ~ richness_std + evenness_std + n_stems_gt5_ha_std + 
    pcoa_1_std + pcoa_2_std + pcoa_3_std + map_std + diurnal_temp_range_std + 
    Matern(1 | x + y), data = dat_std, family = "gaussian")
summary(max_mod)

dd <- dist(dat_std[,c("x","y")])
mm <- MaternCorr(dd, 
  nu = max_mod$corrPars$`1`$nu, rho = max_mod$corrPars$`1`$rho)
plot(as.numeric(dd), as.numeric(mm), 
  xlab = "Distance between pairs of location [in m]", 
  ylab = "Estimated correlation")
sims <- simulateResiduals(max_mod)
plot(sims)

## Predict effect of richness across Zambia 
filled.mapMM(max_mod, add.map = TRUE, 
  plot.axes = quote({axis(1);axis(2)}),
  decorations = quote(points(pred[,coordinates], pch = 17, cex = 0.4)))


dat_std$s1_length_pred <- as.numeric(predict(max_mod, re.form = NA)) 

ggplot() + 
  geom_point(data = dat_std, aes(x = richness_std, y = s1_length)) + 
  geom_smooth(data = dat_std, aes(x = richness_std, y = s1_length_pred))

##' re.form = NA used to remove spatial effects

newdat$calcium <- newdat$calcium + mean(c(0,fixef(m_spamm)[3:4])) # to remove region effect
# get 95% confidence intervals around predictions
newdat <- cbind(newdat, get_intervals(m_spamm, newdata = newdat, intervals = "fixefVar", re.form = NA) + mean(c(0,fixef(m_spamm)[3:4])))


# Define model function
phen_mod <- function(var, pre) {
  # Raw data by group
  pdf(file = paste0("img/z_", pre, "_richness.pdf"), height = 8, width = 10)
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
  message("Model summary")
  print(summary(max_mod))

  # QQ plot
  pdf(file = paste0("img/z_", pre, "_max_mod_qq.pdf"), width = 8, height = 8)
  qqnorm(resid(max_mod))
  qqline(resid(max_mod))
  dev.off()

  # Fixed effect slopes 
  pdf(file = paste0("img/z_", pre, "_max_mod_slopes.pdf"), width = 12, height = 10)
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

  capture.output(mod_compare_anova, file = paste0("out/z_", pre, "_mod_compare.txt"))
}

# Run function for key responses
phen_mod("s1_length", "l")
phen_mod("s1_green_rate", "r")
phen_mod("s1_start", "s")

test_part <- hier.part(dat_std$s1_length, 
  dat_std[,c("richness_std", "evenness_std", "n_stems_gt5_ha_std",
    "pcoa_1_std", "pcoa_2_std", "pcoa_3_std",
    "map_std", "diurnal_temp_range_std")], 
  fam = "gaussian", gof = "Rsqu", barplot = FALSE)

test_rand <- rand.hp(dat_std$s1_length, 
  dat_std[,c("richness_std", "evenness_std", "n_stems_gt5_ha_std",
    "pcoa_1_std", "pcoa_2_std", "pcoa_3_std",
    "map_std", "diurnal_temp_range_std")],
  fam = "gaussian", gof = "Rsqu", num.reps = 10)

# https://www.researchgate.net/publication/282161585_Partitioning_Variance_Into_Constituents_in_Multiple_Regression_Models_Commonality_Analysis

