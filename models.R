# Statistical models 
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-03-01

# Packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggnewscale)
library(shades)
library(xtable)
library(sjPlot)  # devtools::install_github("strengejacke/sjPlot") 
library(sf)
library(MuMIn)
library(ggeffects)
library(emmeans)
library(patchwork)

source("plot_func.R")
source("tex_func.R")

# Import data
plots <- readRDS("dat/plots.rds")
div <- readRDS("dat/div.rds")
modis <- readRDS("dat/modis.rds")
bioclim <- readRDS("dat/bioclim.rds")

# Combine dataframes and scale variables 
dat <- plots %>% 
  inner_join(., div, by = "plot_cluster") %>% 
  inner_join(., modis, by = "plot_cluster") %>% 
  inner_join(., bioclim, by = "plot_cluster") 

greenLagMean <- paste0(round(mean(dat$start_lag, na.rm = TRUE), 0), "$\\pm$", 
  round(sd(dat$start_lag, na.rm = TRUE), 1))

std <- dat %>% 
  mutate_at(
    .vars = names(pred_lookup)[which(names(pred_lookup) != "cluster")],
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  dplyr::select(
    plot_cluster, 
    cluster, 
    ends_with("_std"), 
    names(resp_lookup)) %>%
  rename_at(.vars = vars(ends_with("_std")), 
    .funs = list(~gsub("_std", "", .))) %>% 
  drop_na()

# How many sites after all filtering?
n_sites <- nrow(std)

# MANOVA to show phenological variation within and among vegetation clusters 
dat_manova <- std %>% 
  mutate_at(
    .vars = names(resp_lookup),
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  mutate(cluster = as.character(cluster)) %>% 
  dplyr::select(cluster, ends_with("_std")) %>% 
  st_drop_geometry()

phen_manova <- manova(as.matrix(dat_manova[,grepl("_std", names(dat_manova))]) ~ 
  dat_manova$cluster)

phen_manova_fmt <- paste0("F(", summary(phen_manova)$stats[1], ",", 
  summary(phen_manova)$stats[2], ")=", round(summary(phen_manova)$stats[5], 2), 
  ", ", pFormat(summary(phen_manova)$stats[11]))

# Tukey's tests for each phenological metric per cluster
manova_resp <- names(dat_manova)[grepl("_std", names(dat_manova))]

tukey_out <- do.call(rbind, lapply(manova_resp, function(x) {
  mod <- aov(dat_manova[[x]] ~ dat_manova$cluster)
  mod_means_contr <- emmeans::emmeans(
    object = mod,
    pairwise ~ "cluster",
    adjust = "tukey")
  mod_means <- multcomp::cld(
    object = mod_means_contr$emmeans,
    Letters = letters)
  mod_means$resp <- gsub("_std", "", x)
  return(mod_means)
  }))

tukey_out$.group = trimws(tukey_out$.group)

# Boxplots of phenological metrics per cluster
boxplot_list <- lapply(names(resp_lookup), function(x) {
  tukey_fil <- tukey_out %>% 
    filter(resp == x) %>% 
    mutate(cluster = factor(cluster,
      levels = names(clust_lookup),
      labels = clust_lookup))

  dat_clean <- dat %>%
    st_drop_geometry() %>% 
    mutate(cluster = factor(cluster,
        levels = names(clust_lookup),
        labels = clust_lookup)) 

  ggplot() + 
    geom_boxplot(data = dat_clean,
      aes(x = cluster, y = .data[[x]], fill = cluster)) + 
    scale_fill_manual(name = "Cluster", values = clust_pal) +
    geom_text(data = tukey_fil, 
      aes(x = cluster, 
        y = max(dat_clean[[x]], na.rm = TRUE) + (max(dat_clean[[x]], na.rm = TRUE) * 0.08), 
        label = .group), size = 6) + 
    theme_panel() + 
    theme(axis.text.x = element_blank()) + 
    labs(x = "", y = resp_plot_axes[names(resp_plot_axes) == x])
})

pdf(file = "img/boxplots.pdf", width = 10, height = 8)
wrap_plots(boxplot_list[c(1,3,5,2,4,6)]) +  
  plot_layout(
    ncol = 3, 
    guides = "collect") & theme(legend.position = 'bottom')
dev.off()

# Mixed models of diversity with clusters 

# Define maximal models
precip_vars <- c(
  "cum_precip_seas", 
  "cum_precip_seas", 
  "cum_precip_pre",
  "cum_precip_end", 
  "cum_precip_pre", 
  "cum_precip_seas")

other_vars <- "diurnal_temp_range + eff_rich + diam_quad_mean + Detarioideae + eff_rich:cluster + diam_quad_mean:cluster + cluster"

max_mod_flist <- paste0(names(resp_lookup), " ~ ", precip_vars, " + ", other_vars)

# Fit maximal models
max_ml_list <- lapply(max_mod_flist, function(x) {
  lm(x, data = std, na.action = na.fail)
  })
names(max_ml_list) <- names(resp_lookup)

# Fit models with all combinations of fixed effects and rank by AIC
dredge_list <- lapply(max_ml_list, function(x) {
  dredge(x, evaluate = TRUE, rank = "AIC")
  })
names(dredge_list) <- names(resp_lookup)

# Define climate only "null" models
null_mod_flist <- paste0(names(resp_lookup), " ~ ", 
  precip_vars, " + ", "diurnal_temp_range")

# Fit climate only "null" models
null_ml_list <- lapply(null_mod_flist, function(x) {
  lm(x, data = std)
  })

# Fit REML version of "best" models
# Which models are "best" while including interaction of cluster and richness?
best_mod_flist <- c(
  "EVI_Area ~ cum_precip_seas + diam_quad_mean + Detarioideae + eff_rich + eff_rich:cluster + cluster",
  "season_length ~ cum_precip_seas + eff_rich + diurnal_temp_range + Detarioideae + cluster",
  "green_rate ~ Detarioideae + diurnal_temp_range + eff_rich + cluster",
  "senes_rate ~ Detarioideae + diam_quad_mean + eff_rich + eff_rich:cluster + cluster",
  "start_lag ~ cum_precip_pre + Detarioideae + diurnal_temp_range + eff_rich + eff_rich:cluster + cluster",
  "end_lag ~ diurnal_temp_range + cluster")

best_ml_list <- lapply(best_mod_flist, function(x) {
  lm(x, data = std)
  })

# Get model fit statistics 
modFit <- function(best_ml, null_ml) {
  return(
    data.frame(mod = c("best_ml", "null_ml"),
      aic = c(AIC(best_ml, null_ml)[,2]), 
      bic = c(BIC(best_ml, null_ml)[,2]), 
      rsq = c(summary(best_ml)$adj.r.squared, summary(null_ml)$adj.r.squared),
      logl = c(logLik(best_ml)[1], logLik(null_ml)[1])
    )
  )
}

fit_list <- lapply(seq(length(best_ml_list)), function(x) {
  modFit(best_ml_list[[x]], null_ml_list[[x]])
  })

names(fit_list) <- resp_lookup

# Extract rsq
cum_vi_rsq <- round(fit_list[[1]]$rsq[1], 2)
length_rsq <- round(fit_list[[2]]$rsq[1], 2)
grate_rsq <- round(fit_list[[3]]$rsq[1], 2)
srate_rsq <- round(fit_list[[4]]$rsq[1], 2)
glag_rsq <- round(fit_list[[5]]$rsq[1], 2)
slag_rsq <- round(fit_list[[6]]$rsq[1], 2)

# Nest all model lists in one list
all_mod_list <- list(max_ml_list, dredge_list, null_ml_list, 
  best_ml_list, fit_list)
names(all_mod_list) <- c("max_ml", "dredge", "null_ml",
  "best_ml", "fit_list")

# Model slope plots
mod_slope_df <- do.call(rbind, lapply(seq_along(all_mod_list[[4]]), function(x) {
  mod_pred <- get_model_data(all_mod_list[[4]][[x]], type = "est")

  mod_pred_fil <- mod_pred %>%
    filter(!grepl("cluster", term)) %>%
    mutate(resp = names(resp_lookup)[x]) %>%
    dplyr::select(1:8, group, resp)

  return(mod_pred_fil)
  })) %>%
  mutate(term = factor(term, levels = rev(names(pred_lookup)), labels = rev(pred_lookup)),
    resp = factor(resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)]),
    psig = case_when(
      conf.low >= 0 & conf.high >= 0 | conf.low <= 0 & conf.high <= 0 ~ "*",
      TRUE ~ ""))

pdf(file = "img/mod_slopes.pdf", width = 10, height = 5)
ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(data = mod_slope_df, 
    aes(xmin = conf.low, xmax = conf.high, y = term, colour = group),
    height = 0) + 
  geom_point(data = mod_slope_df,
    aes(x = estimate, y = term, fill = group),
    shape = 21, colour = "black") + 
  geom_text(data = mod_slope_df,
    aes(x = estimate, y = term, colour = group, label = psig),
    size = 8, nudge_y = 0.1) + 
  facet_wrap(~resp, scales = "free_x") + 
  scale_fill_manual(name = "", values = pal[3:4]) +
  scale_colour_manual(name = "", values = pal[3:4]) +
  theme_panel() + 
  theme(legend.position = "none") + 
  labs(x = "Standardised slope coefficient", y = "")
dev.off()

# Export interaction effect on richness plots
intf <- function(x) {
  best_int_mod_flist <- best_mod_flist[grepl(x, best_mod_flist)] %>%
    ifelse(grepl(paste0(x, ":cluster"), .), ., paste0(., " + ", x, ":cluster"))

  best_int_ml_list <- lapply(best_int_mod_flist, function(x) {
    lm(formula(x), data = std)
    })

  marg_df <- do.call(rbind, lapply(seq_along(best_int_ml_list), function(y) {
    preds <- ggemmeans(best_int_ml_list[[y]], terms = c(x, "cluster"), 
      interval = "confidence", ci.lvl = 0.95)
    preds$resp <- gsub("\\s~.*", "", best_int_mod_flist[[y]])
    preds$pred_name <- x
    names(preds) <- c("val", "pred", "se", "conf_lo", "conf_hi", "cluster", "resp", "pred_name")
    return(preds)
  }))

  return(list(best_int_ml_list, marg_df))
}

eff_rich_intf <- intf("eff_rich")
diam_quad_mean_intf <- intf("diam_quad_mean")
Detariodeae_intf <- intf("Detarioideae")

marg_df <- rbind(eff_rich_intf[[2]], diam_quad_mean_intf[[2]], Detariodeae_intf[[2]])

marg_df$resp_clean <- factor(marg_df$resp, 
      levels = names(resp_plot_axes[c(1,3,5,2,4,6)]), 
      labels = resp_plot_axes[c(1,3,5,2,4,6)])

marg_df$pred_name_clean <- factor(marg_df$pred_name, 
      levels = names(pred_lookup), 
      labels = pred_lookup)

marg_df$cluster_clean <- factor(marg_df$cluster,
      levels = names(clust_lookup),
      labels = clust_lookup)

pdf(file = "img/mod_marg.pdf", width = 10, height = 11)
ggplot() + 
  geom_ribbon(data = marg_df, 
    aes(x = val, ymin = conf_lo, ymax = conf_hi, colour = cluster_clean), 
    linetype = 2, alpha = 0.2, fill = NA) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) +
  new_scale_colour() + 
  geom_line(data = marg_df, 
    aes(x = val, y = pred, colour = cluster_clean),
    linewidth = 1.5) + 
  scale_colour_manual(name = "Cluster", 
    values = brightness(clust_pal, 0.75)) +
  facet_grid(resp_clean~pred_name_clean, scales = "free") + 
  theme_panel() + 
  theme(legend.position = "bottom") + 
  labs(x = "", y = "")
dev.off()

# Write variables
write(
  c(
    commandOutput(phen_manova_fmt, "phenManova"),
    commandOutput(n_sites, "nSites"),
		commandOutput(cum_vi_rsq, "cumviRsq"),
		commandOutput(length_rsq, "lengthRsq"),
		commandOutput(grate_rsq, "grateRsq"),
		commandOutput(srate_rsq, "srateRsq"),
		commandOutput(glag_rsq, "glagRsq"),
		commandOutput(slag_rsq, "slagRsq"),
    commandOutput(greenLagMean, "greenLagMean")
    ),
  file = "out/models_vars.tex")

