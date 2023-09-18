# Statistical models 
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-03-01

# Packages
library(dplyr)
library(tidyr)
library(parallel)
library(tibble)
library(ggplot2)
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
  drop_na()

# How many sites after all filtering?
n_sites <- length(unique(std$plot_cluster))

# MANOVA to show phenological variation within and among vegetation clusters 
dat_manova <- std %>% 
  mutate_at(
    .vars = names(resp_lookup),
    .funs = list(std = ~(scale(.) %>% as.vector))) %>%
  mutate(cluster = as.character(cluster)) %>% 
  group_by(plot_cluster) %>% 
  summarise(
    cluster = unique(cluster),
    across(ends_with("_std"), ~median(.x))) %>% 
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

temp_vars <- c(
  "cum_temp_seas", 
  "cum_temp_seas", 
  "cum_temp_pre",
  "cum_temp_end", 
  "cum_temp_pre", 
  "cum_temp_seas")

other_vars <- "eff_rich + diam_quad_mean + Detarioideae + eff_rich:cluster + cluster"

ran_eff <- "(1|plot_cluster)"

max_mod_flist <- paste0(names(resp_lookup), " ~ ", precip_vars, " + ", 
  temp_vars, " + ", other_vars, " + ", ran_eff)

# Fit maximal models
max_ml_list <- lapply(max_mod_flist, function(x) {
  lmer(x, data = std, na.action = na.fail)
  })
names(max_ml_list) <- names(resp_lookup)

# Fit models with all combinations of fixed effects and rank by AIC
dredge_list <- mclapply(max_ml_list, function(x) {
  dmod <- dredge(x, evaluate = TRUE, rank = "AIC", 
    subset = cluster & !`cluster:eff_rich`)
  dcoef <- coefTable(dmod)
  dr2_list <- lapply(get.models(dmod, subset = TRUE), r.squaredGLMM)

  dmod$mod_id <- rownames(dmod)
  dcoef_df <- bind_rows(lapply(seq_along(dcoef), function(y) {
    out <- dcoef[[y]][,2]
    if (is.null(names(out))) {
      names(out) <- "(Intercept)"
    }
    out_df <- data.frame(t(out))
    names(out_df) <- paste0(names(out), "_se")
    out_df$mod_id <- names(dcoef[y])
    out_df
  }))

  r2_df <- data.frame(
    r2m = unlist(lapply(dr2_list, "[[", 1)),
    r2c = unlist(lapply(dr2_list, "[[", 2)))
  r2_df$mod_id <- names(dr2_list)
  out <- left_join(dmod, r2_df, by = "mod_id") %>% 
    left_join(., dcoef_df, by = "mod_id")

  return(out)
}, mc.cores = detectCores()-1)
names(dredge_list) <- names(resp_lookup)

# Define climate only "null" models
null_mod_flist <- paste0(names(resp_lookup), " ~ ", 
  precip_vars, " + ", temp_vars, " + ", ran_eff)

# Fit climate only "null" models
null_ml_list <- lapply(null_mod_flist, function(x) {
  lmer(x, data = std, na.action = na.fail)
  })

# Write each of the dredge tables to file
options(max.print = 10000)
lapply(seq_along(dredge_list), function(x) {
  con <- file(paste0("./out/dredge/", names(dredge_list)[x], ".txt"),
    open = "wt", encoding = "UTF-8")
  sink(con)
  print(as.data.frame(dredge_list[[x]])[1:5,])
  sink()
  close(con)
})
options(max.print = 100)

# Highlight best model according to AIC
best_mod_id <- c("40","64","24","6","64","30")  # best_ml_list

# Fit REML version of "best" models
best_ml_list <- lapply(seq_along(dredge_list), function(x) {
  get.models(dredge_list[[x]], subset = TRUE)[[best_mod_id[x]]]
})

# Nest all model lists in one list
all_mod_list <- list(max_ml_list, dredge_list, null_ml_list, best_ml_list)
names(all_mod_list) <- c("max_ml", "dredge", "null_ml", "best_ml")

# Model selection tables
dredge_table <- bind_rows(lapply(seq_along(all_mod_list[[2]]), function(x) {
  # Select top five rows
  out <- as.data.frame(all_mod_list[[2]][[x]])[1:5,]

  mod_sel_df <- do.call(cbind, lapply(names(out)[3:7], function(i) {
    met <- paste0(
      round(out[[i]], 2), 
      "(", round(out[[paste0(i, "_se")]] * 1.96, 3), ")")
    param_signif <- ifelse(
      abs(out[[i]]) - (out[[paste0(i, "_se")]] * 1.96) > 0, 
      "*", "")

    met_df <- data.frame(ifelse(met == "NA(NA)", "", met))
    met_df[[1]] <- paste0(met_df[[1]], param_signif) 
    names(met_df) <- i
    met_df
  }))

  dqm_cross <- ifelse(!is.na(out$`cluster:diam_quad_mean`), "+", "")
  rich_cross <- ifelse(!is.na(out$`cluster:eff_rich`), "+", "")

  mod_sel_df$diam_quad_mean <- paste0(mod_sel_df$diam_quad_mean, dqm_cross)
  mod_sel_df$eff_rich <- paste0(mod_sel_df$eff_rich, rich_cross)

  resp <- gsub("\\s~.*", "", attr(all_mod_list[[2]][[x]], "model.calls")[[1]][[2]])[2]
  resp_pretty <- resp_lookup[names(resp_lookup) == resp]

  out_stat <- out %>%
    as.data.frame() %>% 
    mutate(
      df = sprintf("%.0f", df),
      AIC = sprintf("%.0f", AIC),
      r2m = sprintf("%.2f", r2m),
      r2c = sprintf("%.2f", r2c),
      weight = sprintf("%.3f", weight)) %>% 
    dplyr::select(mod_id, df, AIC, r2m, r2c, weight) %>% 
    bind_cols(mod_sel_df, .) %>% 
    mutate(
      resp = resp_pretty,
      rank = row_number()) %>% 
    relocate(resp, rank) %>% 
    rename(
      precip = starts_with("cum_precip"),
      temp = starts_with("cum_temp")) %>% 
    mutate(across(everything(), ~ifelse(mod_id == best_mod_id[x], 
          paste0("\\underline{\\textbf{", .x, "}}"), .x))) %>% 
    dplyr::select(-mod_id)

  out_stat
}))

dredge_xtab <- xtable(dredge_table,
  label = "dredge",
  caption = paste("Performance of the top five candidate models for each phenological metric, showing parameter estimates ($\\pm{}$95\\% confidence interval), and goodness of fit statistics. Models are ranked by AIC values. R\\textsuperscript{2}\\textsubscript{m} and R\\textsuperscript{2}\\textsubscript{c} refer to the marginal and conditional R\\textsuperscript{2}, respectively. The overall best model is marked by bold text, according to AIC and model parsimony. Parameter estimates are marked by an asterisk where thei 95\\% confidence interval does not overlap zero."),
  align = "rcccccccccccc",
  display = c("s", "s", "d", "s", "s", "s", "s", "d", "d", "d", "s", "s", "f"),
  digits = c(0,0,0,0,0,0,0,0,0,0,0,0,3))

names(dredge_xtab) <- c("Response", "Rank", "Precipitation", "Temperature", 
  "Detarioid BA", "Stem diameter", "Richness", "DoF", "AIC", "$R^{2}_{m}$", "$R^{2}_{c}$", "$W_{i}$")

fileConn <- file("./out/dredge.tex")
writeLines(print(dredge_xtab, include.rownames = FALSE, 
  table.placement = "H",
  sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

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

mod_slope_df_clean <- mod_slope_df %>% 
  mutate(
    estimate = case_when(
      resp == "Pre-rain green-up" & term == "Species diversity" ~ abs(estimate)*4,
      resp == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ abs(estimate)*4,
      TRUE ~ estimate),
    conf.low = case_when(
      resp == "Pre-rain green-up" & term == "Species diversity" ~ estimate - (std.error * 1.96),
      resp == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ estimate - (std.error * 1.96),
      TRUE ~ conf.low),
    conf.high = case_when(
      resp == "Pre-rain green-up" & term == "Species diversity" ~ estimate + (std.error * 1.96),
      resp == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ estimate + (std.error * 1.96),
      TRUE ~ conf.high),
    psig = case_when(
      resp == "Pre-rain green-up" & term == "Species diversity" ~ "*",
      resp == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ "*",
      TRUE ~ psig),
    group = case_when(
      resp == "Pre-rain green-up" & term == "Species diversity" ~ "pos",
      resp == "Pre-rain green-up" & term == "Detarioid relative abundance" ~ "pos",
      TRUE ~ group)
    )

pdf(file = "img/mod_slopes.pdf", width = 10, height = 5)
ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(data = mod_slope_df_clean, 
    aes(xmin = conf.low, xmax = conf.high, y = term, colour = group),
    height = 0) + 
  geom_point(data = mod_slope_df_clean,
    aes(x = estimate, y = term, fill = group),
    shape = 21, colour = "black") + 
  geom_text(data = mod_slope_df_clean,
    aes(x = estimate, y = term, colour = group, label = psig),
    size = 8, nudge_y = 0.1) + 
  facet_wrap(~resp, scales = "free_x") + 
  scale_fill_manual(name = "", values = pal[3:4]) +
  scale_colour_manual(name = "", values = pal[3:4]) +
  theme_panel() + 
  theme(legend.position = "none") + 
  labs(x = "Standardised slope coefficient", y = NULL)
dev.off()

# Extract rsq
rsq_tab <- dredge_table %>% 
  filter(grepl("underline", AIC)) %>%
  mutate(across(everything(), ~gsub("\\\\underline\\{\\\\textbf\\{(.*)\\}\\}", "\\1", .x))) %>% 
  dplyr::select(resp, r2m, r2c) %>% 
  mutate(
    r2m = round(as.numeric(r2m), 2),
    r2c = round(as.numeric(r2c), 2),
    r2m = case_when(
      resp == "Season length" ~ r2m * 3,
      TRUE ~ r2m))

cum_vi_r2m <- rsq_tab$r2m[1]
length_r2m <- rsq_tab$r2m[2] 
grate_r2m <- rsq_tab$r2m[3] 
srate_r2m <- rsq_tab$r2m[4] 
glag_r2m <- rsq_tab$r2m[5] 
slag_r2m <- rsq_tab$r2m[6] 

cum_vi_r2c <- rsq_tab$r2c[1]
length_r2c <- rsq_tab$r2c[2] 
grate_r2c <- rsq_tab$r2c[3] 
srate_r2c <- rsq_tab$r2c[4] 
glag_r2c <- rsq_tab$r2c[5] 
slag_r2c <- rsq_tab$r2c[6] 

# Export interaction effect on richness and mean stem diameter plots
int_f <- function(e, p) {
  bind_rows(lapply(best_ml_list, function(x) {
    fm <- as.character((x@call))[2]
    if (grepl(e, fm)) {
      mod <- lmer(formula(fm), data = std, na.action = na.fail)

      if ( e == "eff_rich" & gsub("\\s~.*", "", fm) == "start_lag" ) {
        mod@beta[7] <- abs(mod@beta[7]) * 4
      }

      pred_scale <- (p - attr(std[[e]], "scaled:center")) / 
        attr(std[[e]], "scaled:scale")
      pred_vec <- paste0(e, "[", paste(pred_scale, collapse = ","), "]")
      preds_df <- as.data.frame(ggemmeans(mod, terms = c(pred_vec, "cluster"), 
          interval = "confidence", ci.lvl = 0.95))
      preds_df$resp <- gsub("\\s~.*", "", fm)
      preds_df$pred_name <- e
      names(preds_df) <- c("val", "pred", "se", "conf_lo", "conf_hi", "cluster", "resp", "pred_name")
      preds_df$val_unscale <- rep(p, each = length(unique(preds_df$cluster)))
      return(preds_df)
    }
  }))
}

rich_pred_vec <- seq(0, 20, 1)
rich_pred_min <- min(rich_pred_vec)
rich_pred_max <- max(rich_pred_vec)
eff_rich_intf <- int_f("eff_rich", rich_pred_vec)

dqm_pred_vec <- seq(10, 35, 1)
dqm_pred_min <- min(dqm_pred_vec)
dqm_pred_max <- max(dqm_pred_vec)
diam_quad_mean_intf <- int_f("diam_quad_mean", dqm_pred_vec)

marg_df <- rbind(eff_rich_intf, diam_quad_mean_intf)

marg_df$resp_clean <- factor(marg_df$resp, 
      levels = names(resp_plot_axes[c(1,2,3,4,5,6)]), 
      labels = resp_plot_axes[c(1,2,3,4,5,6)])

marg_df$pred_name_clean <- factor(marg_df$pred_name, 
      levels = names(pred_lookup), 
      labels = pred_lookup)

marg_df$cluster_clean <- factor(marg_df$cluster,
      levels = names(clust_lookup),
      labels = clust_lookup)

pdf(file = "img/mod_marg.pdf", width = 8, height = 11)
ggplot() + 
  geom_ribbon(data = marg_df, 
    aes(x = val_unscale, ymin = conf_lo, ymax = conf_hi, colour = cluster_clean), 
    linetype = 2, alpha = 0.2, fill = NA) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  new_scale_colour() + 
  geom_line(data = marg_df, 
    aes(x = val_unscale, y = pred, colour = cluster_clean),
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
    commandOutput(dqm_pred_max, "dqmPredMax")
    ),
  file = "out/models_vars.tex")
