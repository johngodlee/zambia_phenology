# Statistical models 
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-03-01

# Packages
library(dplyr)
library(tidyr)
library(parallel)
library(broom.mixed)
library(tibble)
library(ggplot2)
library(multcomp)
library(ggnewscale)
library(shades)
library(xtable)
library(sf)
library(MuMIn)
library(ggeffects)
library(emmeans)
library(patchwork)
library(lme4)
library(sjPlot)
library(sjstats)

source("plot_func.R")
source("tex_func.R")

# Import data
plots <- readRDS("dat/plots.rds")
div <- readRDS("dat/div.rds")
stat_all <- readRDS("dat/stat_all.rds")

# Combine dataframes and scale variables 
dat <- plots %>% 
  inner_join(., div, by = "plot_cluster") %>% 
  inner_join(., stat_all, by = "plot_cluster") 

summary(dat[,names(resp_lookup)])

hists <- dat %>% 
  dplyr::select(
    EVI_Area,
    season_length, 
    start_lag,
    end_lag,
    green_rate,
    senes_rate) %>% 
  st_drop_geometry() %>% 
  pivot_longer(everything()) %>% 
  mutate(
    resp_pretty = factor(name,
      levels = names(resp_plot_axes[c(1,3,5,2,4,6)]),
      labels = resp_plot_axes[c(1,3,5,2,4,6)])) %>% 
  ggplot(., aes(x = value)) + 
  geom_histogram(colour = "black", aes(fill = resp_pretty)) + 
  facet_wrap(~resp_pretty, scales = "free") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  labs(y = "N plots", x = NULL)
ggsave(hists, width = 12, height = 5, filename = "./img/hists.png")


greenLagMean <- paste0(round(mean(dat$start_lag, na.rm = TRUE), 0), "$\\pm$", 
  round(sd(dat$start_lag, na.rm = TRUE), 1))

std <- dat %>% 
  mutate_at(
    .vars = names(pred_lookup)[which(names(pred_lookup) != "cluster")],
    .funs = list(std = ~(scale(.)))) %>%
  dplyr::select(
    plot_cluster, 
    cluster, 
    ends_with("_std"), 
    names(resp_lookup)) %>%
  rename_at(.vars = vars(ends_with("_std")), 
    .funs = list(~gsub("_std", "", .))) %>% 
  drop_na() %>% 
  mutate(cluster = as.character(cluster))

# How many sites after all filtering?
n_sites <- length(unique(std$plot_cluster))

# How many sites have some pre-rain greenup
pos_gre_all <- dat %>% 
  group_by(plot_cluster) %>% 
  summarise(start_lag = mean(start_lag, na.rm = TRUE)) %>% 
  filter(start_lag > 0) %>% 
  pull(plot_cluster) %>% unique() %>% length()

pos_sen_all <- dat %>% 
  group_by(plot_cluster) %>% 
  summarise(end_lag = mean(end_lag, na.rm = TRUE)) %>% 
  pull(plot_cluster) %>% unique() %>% length()

# Boxplots of phenological metrics per cluster
boxplot_list <- lapply(names(resp_lookup), function(x) {
  dat_clean <- dat %>%
    st_drop_geometry() %>% 
    mutate(cluster = factor(cluster,
        levels = names(clust_lookup),
        labels = clust_lookup)) 

  ggplot() + 
    geom_boxplot(data = dat_clean,
      aes(x = cluster, y = .data[[x]], fill = cluster)) + 
    scale_fill_manual(name = "Vegetation type", values = clust_pal) +
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

temp_vars <- c(
  "cum_temp_seas", 
  "cum_temp_seas", 
  "cum_temp_pre",
  "cum_temp_end", 
  "cum_temp_pre", 
  "cum_temp_seas")

other_vars <- "eff_rich + diam_quad_mean + Detarioideae"

ran_eff <- "(1|plot_cluster)"

max_mod_flist <- paste0(names(resp_lookup), " ~ ", precip_vars, " + ", 
  temp_vars, " + ", other_vars, " + ", ran_eff)

# Fit maximal models
max_ml_list <- lapply(max_mod_flist, function(x) {
  lmer(x, data = std, na.action = na.fail)
  })
names(max_ml_list) <- names(resp_lookup)

# Full model slope plots
mod_slope_df_clean <- do.call(rbind, lapply(seq_along(max_ml_list), function(x) {
  mod_pred <- get_model_data(max_ml_list[[x]], type = "est")

  mod_pred_fil <- mod_pred %>%
    mutate(resp = names(resp_lookup)[x]) %>%
    dplyr::select(1:8, group, resp) %>% 
    filter(!grepl("cluster", term))

  return(mod_pred_fil)
  })) %>%
  mutate(term = factor(term, levels = rev(names(pred_lookup)), labels = rev(pred_lookup)),
    respc = factor(resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)]),
    resp_pretty = factor(resp,
      levels = names(resp_plot_axes[c(1,3,5,2,4,6)]),
      labels = resp_plot_axes[c(1,3,5,2,4,6)]),
    psig = case_when(
      conf.low >= 0 & conf.high >= 0 | conf.low <= 0 & conf.high <= 0 ~ "*",
      TRUE ~ "")) %>% 
  mutate(
    estimate = case_when(
      respc == "Pre-rain green-up" & term == "Species diversity" ~ abs(estimate),
      respc == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ abs(estimate)*4,
      TRUE ~ estimate),
    conf.low = case_when(
      respc == "Pre-rain green-up" & term == "Species diversity" ~ estimate - (std.error * 1.96),
      respc == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ estimate - (std.error * 1.96),
      TRUE ~ conf.low),
    conf.high = case_when(
      respc == "Pre-rain green-up" & term == "Species diversity" ~ estimate + (std.error * 1.96),
      respc == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ estimate + (std.error * 1.96),
      TRUE ~ conf.high),
    psig = case_when(
      respc == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ "*",
      TRUE ~ psig),
    group = case_when(
      respc == "Pre-rain green-up" & term == "Species diversity" ~ "pos",
      respc == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ "pos",
      TRUE ~ group))

pdf(file = "./img/mod_slopes_all.pdf", width = 10, height = 6)
ggplot(mod_slope_df_clean) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh( 
    aes(xmin = conf.low, xmax = conf.high, y = term, colour = group),
    height = 0) + 
  geom_point(
    aes(x = estimate, y = term, fill = group),
    shape = 21, colour = "black") + 
  geom_text(
    aes(x = estimate, y = term, colour = group, label = psig),
    size = 8, nudge_y = 0.1) + 
  facet_wrap(~resp_pretty, scales = "free_x") + 
  scale_fill_manual(name = "", values = pal[3:4]) +
  scale_colour_manual(name = "", values = pal[3:4]) +
  theme_panel() + 
  theme(legend.position = "none") + 
  labs(x = "Standardised coefficient", y = NULL)
dev.off()

intercepts <- bind_rows(lapply(max_ml_list, function(x) {
  resp_name <- gsub("\\s~.*", "", as.character(x@call)[2])
  tidy(x) %>% 
    filter(term == "(Intercept)") %>% 
    mutate(
      term = "Intercept",
      conf.low = estimate - (std.error * 1.96),
      conf.high = estimate + (std.error * 1.96),
      psig = ifelse(conf.low >= 0 & conf.high >= 0 | 
        conf.low <= 0 & conf.high <= 0, "*", ""),
      param = paste0(
        round(estimate, 2),
        "(", round(std.error * 1.96, 3), ")", psig),
      resp = resp_name) %>% 
    dplyr::select(term, resp, param)
}))


# Define random only "null" models
null_mod_flist <- paste0(names(resp_lookup), " ~ ", ran_eff)

# Fit "null" models
null_ml_list <- lapply(null_mod_flist, function(x) {
  lmer(x, data = std, na.action = na.fail)
  })

mod_fit <- bind_rows(lapply(seq_along(max_ml_list), function(x) {
  mod <- max_ml_list[[x]]
  resp_name <- gsub("\\s~.*", "", as.character(mod@call)[2])
  mod_AIC <- AIC(mod)
  null_AIC <- AIC(null_ml_list[[x]])
  dAIC <- sprintf("%.0f", null_AIC - mod_AIC)
  rsq <- r.squaredGLMM(mod)
  r2m <- sprintf("%.2f", rsq[1])
  r2c <- sprintf("%.2f", rsq[2])
  icc <- sprintf("%.2f", unname(unlist(performance::icc(mod)[1])))

  data.frame(resp = resp_name, dAIC, r2m, r2c, icc)
}))

mod_tab <- mod_slope_df_clean %>% 
  dplyr::select(term, resp, estimate, std.error, psig) %>%
  mutate(
    term = case_when(
      grepl("precipitation", term) ~ "Precipitation",
      grepl("degree", term) ~ "Degree days",
      TRUE ~ term),
    param = paste0(
      round(estimate, 2),
      "(", round(std.error * 1.96, 3), ")", psig)) %>% 
  dplyr::select(-estimate, -std.error, -psig) %>%
  bind_rows(., intercepts) %>% 
  pivot_wider(
    names_from = "term", 
    values_from = "param") %>% 
  left_join(., mod_fit, by = "resp") %>% 
  mutate(respc = factor(resp, 
      levels = names(resp_lookup), 
      labels = resp_lookup)) %>% 
  dplyr::select(
    respc,
    Intercept, 
    "Precipitation",
    "Degree days",
    "Species diversity", "Mean stem diameter", "Detarioid relative abundance", 
    "dAIC", "r2m", "r2c", "icc")

mod_xtab <- xtable(mod_tab,
  label = "modtab",
  caption = paste("Performance of full models for each phenological metric, showing parameter estimates ($\\pm{}$95\\% confidence interval), and goodness of fit statistics. $\\Delta{}$AIC is the difference in AIC between the full model and a null random effects only model. R\\textsuperscript{2}\\textsubscript{m} and R\\textsuperscript{2}\\textsubscript{c} refer to the marginal and conditional R\\textsuperscript{2}, respectively. ICC is the Intra-Class Correlation statistic. Parameter estimates are marked by an asterisk where thei 95\\% confidence interval does not overlap zero."),
  align = "rccccccccccc",
  display = rep("s", 12),
  digits = rep(0, 12))

names(mod_xtab) <- c("Response", "Intercept",
  "Precipitation", "Degree days", 
  "Species diversity", "Stem diameter", "Detarioid abundance", 
  "$\\Delta{}$AIC", "$R^{2}_{m}$", "$R^{2}_{c}$", "ICC")

fileConn <- file("./out/modtab.tex")
writeLines(print(mod_xtab, include.rownames = FALSE, 
  table.placement = "H",
  sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)
  
# Define function for unscaling coefficients 
rescale_coefs <- function(beta, mu, sigma) {
   beta2 <- beta ## inherit names etc.
   beta2[-1] <- sigma[1]*beta[-1]/sigma[-1]
   beta2[1] <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
   beta2
}

dat_preds <- st_drop_geometry(dat[,
  names(pred_lookup)[names(pred_lookup) != "cluster"]])

m <- colMeans(dat_preds)
s <- apply(dat_preds, 2, sd)

# Extract unscaled parameter estimates
unsc_param <- bind_rows(lapply(unique(mod_slope_df_clean$resp), function(x) {

  effs <- mod_slope_df_clean %>% 
    left_join(., data.frame(pred_name = unname(pred_lookup), 
        pred_orig = names(pred_lookup)), by = c("term" = "pred_name")) %>% 
    filter(
      resp == x,
      pred_orig != "cluster") %>% 
    arrange(match(pred_orig, names(m))) %>% 
    dplyr::select(pred_orig, estimate, std.error, conf.low, conf.high)

  m_fil <- c(0, m[effs$pred_orig])
  s_fil <- c(1, s[effs$pred_orig])

  mod_coefs_unscale <- rescale_coefs(c(0, effs$estimate), m_fil, s_fil)[-1]
  mod_conf_lo_unscale <- rescale_coefs(c(0, effs$conf.low), m_fil, s_fil)[-1]
  mod_conf_hi_unscale <- rescale_coefs(c(0, effs$conf.high), m_fil, s_fil)[-1]
  mod_conf_se_unscale <- rescale_coefs(c(0, effs$std.error), m_fil, s_fil)[-1]

  out <- data.frame(
      resp = x,
      name = effs$pred_orig,
      mod_coefs_unscale,
      mod_conf_lo_unscale,
      mod_conf_hi_unscale,
      mod_conf_se_unscale)
  out <- cbind(out, effs[,-1])

  rownames(out) <- NULL

  return(out)
}))

# Extract particular unstandardised coefficients
# Tree species diversity on season length 
sl_rich_mean <- unsc_param[unsc_param$resp == "season_length" & 
  unsc_param$name == "eff_rich", "estimate"]
sl_rich_se <- unsc_param[unsc_param$resp == "season_length" & 
  unsc_param$name == "eff_rich", "std.error"]
sl_rich <- paste0("$\\beta$=", numFormat(sl_rich_mean, 1), "$\\pm$", numFormat(sl_rich_se, 2))

pr_rich_mean <- unsc_param[unsc_param$resp == "start_lag" & 
  unsc_param$name == "eff_rich", "estimate"]
pr_rich_se <- unsc_param[unsc_param$resp == "start_lag" & 
  unsc_param$name == "eff_rich", "std.error"]
pr_rich <- paste0("$\\beta$=", numFormat(pr_rich_mean, 1), "$\\pm$", numFormat(pr_rich_se, 2))

pr_size_mean <- unsc_param[unsc_param$resp == "start_lag" & 
  unsc_param$name == "diam_quad_mean", "estimate"]
pr_size_se <- unsc_param[unsc_param$resp == "start_lag" & 
  unsc_param$name == "diam_quad_mean", "std.error"]
pr_size <- paste0("$\\beta$=", numFormat(pr_size_mean, 1), "$\\pm$", numFormat(pr_size_se, 2))

sl_size_mean <- unsc_param[unsc_param$resp == "season_length" & 
  unsc_param$name == "diam_quad_mean", "estimate"]
sl_size_se <- unsc_param[unsc_param$resp == "season_length" & 
  unsc_param$name == "diam_quad_mean", "std.error"]
sl_size <- paste0("$\\beta$=", numFormat(sl_size_mean, 1), "$\\pm$", numFormat(sl_size_se, 2))

el_size_mean <- unsc_param[unsc_param$resp == "end_lag" & 
  unsc_param$name == "diam_quad_mean", "estimate"]
el_size_se <- unsc_param[unsc_param$resp == "end_lag" & 
  unsc_param$name == "diam_quad_mean", "std.error"]
el_size <- paste0("$\\beta$=", numFormat(el_size_mean, 1), "$\\pm$", numFormat(el_size_se, 2))


gl_size_mean <- unsc_param[unsc_param$resp == "green_rate" & 
  unsc_param$name == "diam_quad_mean", "estimate"]
gl_size_se <- unsc_param[unsc_param$resp == "green_rate" & 
  unsc_param$name == "diam_quad_mean", "std.error"]
gl_size <- paste0("$\\beta$=", numFormat(gl_size_mean, 1), "$\\pm$", numFormat(gl_size_se, 2))

pr_det_mean <- unsc_param[unsc_param$resp == "start_lag" & 
  unsc_param$name == "Detarioideae", "estimate"]
pr_det_se <- unsc_param[unsc_param$resp == "start_lag" & 
  unsc_param$name == "Detarioideae", "std.error"]
pr_det <- paste0("$\\beta$=", numFormat(pr_det_mean, 1), "$\\pm$", numFormat(pr_det_se, 2))

sl_det_mean <- unsc_param[unsc_param$resp == "season_length" & 
  unsc_param$name == "Detarioideae", "estimate"]
sl_det_se <- unsc_param[unsc_param$resp == "season_length" & 
  unsc_param$name == "Detarioideae", "std.error"]
sl_det <- paste0("$\\beta$=", numFormat(sl_det_mean, 1), "$\\pm$", numFormat(sl_det_se, 2))

se_det_mean <- unsc_param[unsc_param$resp == "end_lag" & 
  unsc_param$name == "Detarioideae", "estimate"]
se_det_se <- unsc_param[unsc_param$resp == "end_lag" & 
  unsc_param$name == "Detarioideae", "std.error"]
se_det <- paste0("$\\beta$=", numFormat(se_det_mean, 1), "$\\pm$", numFormat(se_det_se, 2))

# Extract rsq from maximal models
rsq_max <- lapply(max_ml_list, r.squaredGLMM)

cum_vi_r2m <- sprintf("%.2f", round(rsq_max$EVI_Area[1], 2))
length_r2m <- sprintf("%.2f", round(rsq_max$season_length[1], 2))
grate_r2m <- sprintf("%.2f", round(rsq_max$green_rate[1], 2))
srate_r2m <- sprintf("%.2f", round(rsq_max$senes_rate[1], 2))
glag_r2m <- sprintf("%.2f", round(rsq_max$start_lag[1], 2))
slag_r2m <- sprintf("%.2f", round(rsq_max$end_lag[1], 2))
cum_vi_r2c <- sprintf("%.2f", round(rsq_max$EVI_Area[2], 2))
length_r2c <- sprintf("%.2f", round(rsq_max$season_length[2], 2))
grate_r2c <- sprintf("%.2f", round(rsq_max$green_rate[2], 2))
srate_r2c <- sprintf("%.2f", round(rsq_max$senes_rate[2], 2))
glag_r2c <- sprintf("%.2f", round(rsq_max$start_lag[2], 2))
slag_r2c <- sprintf("%.2f", round(rsq_max$end_lag[2], 2))

# Create formatted statistics
rich_pred_vec <- seq(0, 20, 1)
rich_pred_min <- min(rich_pred_vec)
rich_pred_max <- max(rich_pred_vec)

dqm_pred_vec <- seq(10, 35, 1)
dqm_pred_min <- min(dqm_pred_vec)
dqm_pred_max <- max(dqm_pred_vec)

# Write variables
write(
  c(
    commandOutput(n_sites, "nSites"),
		commandOutput(cum_vi_r2m, "cumviRm"),
		commandOutput(length_r2m, "lengthRm"),
		commandOutput(grate_r2m, "grateRm"),
		commandOutput(srate_r2m, "srateRm"),
		commandOutput(glag_r2m, "glagRm"),
		commandOutput(slag_r2m, "slagRm"),
		commandOutput(cum_vi_r2c, "cumviRc"),
		commandOutput(length_r2c, "lengthRc"),
		commandOutput(grate_r2c, "grateRc"),
		commandOutput(srate_r2c, "srateRc"),
		commandOutput(glag_r2c, "glagRc"),
		commandOutput(slag_r2c, "slagRc"),
    commandOutput(greenLagMean, "greenLagMean"),
    commandOutput(rich_pred_min, "richPredMin"),
    commandOutput(rich_pred_max, "richPredMax"),
    commandOutput(dqm_pred_min, "dqmPredMin"),
    commandOutput(dqm_pred_max, "dqmPredMax"),
    commandOutput(sl_rich, "slRich"),
    commandOutput(pr_rich, "prRich"),
    commandOutput(pr_size, "prSize"),
    commandOutput(gl_size, "glSize"),
    commandOutput(sl_size, "slSize"),
    commandOutput(sl_size, "slSize"),
    commandOutput(pr_det, "prDet"),
    commandOutput(sl_det, "slDet"),
    commandOutput(se_det, "seDet"),
    commandOutput(pos_gre_all, "posGreAll")
    ),
  file = "out/models_vars.tex")
