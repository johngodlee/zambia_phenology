# Testing effect of species composition and richness on land surface phenology with Zambia ILUAii data 
# John Godlee (johngodlee@gmail.com)
# 2020-10-25

# Set working directory

# Packages
library(dplyr)
library(nlme)
library(MuMIn)
library(sjPlot)
library(ggeffects)
library(ggplot2)
library(ggnewscale)
library(shades)
library(xtable)
library(car)

source("functions.R")

# Import data
dat_std <- readRDS("dat/plots_anal.rds")

# Mixed models with cluster random effects ----

# Fit candidate models
max_ml_c <- gls(cum_vi ~ cum_precip_seas + diurnal_temp_range + 
  evenness + richness + richness:cluster, 
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")
dredge_c <- dredge(max_ml_c, evaluate = TRUE, rank = "AIC")

max_ml_l <- gls(s1_length ~ cum_precip_seas + diurnal_temp_range + 
  evenness + richness + richness:cluster,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")
dredge_l <- dredge(max_ml_l, evaluate = TRUE, rank = "AIC")

max_ml_r <- gls(s1_green_rate ~ cum_precip_pre + diurnal_temp_range +
  evenness + richness + richness:cluster,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")
dredge_r <- dredge(max_ml_r, evaluate = TRUE, rank = "AIC")

max_ml_d <- gls(s1_senes_rate ~ cum_precip_end + diurnal_temp_range + 
  evenness + richness + richness:cluster,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")
dredge_d <- dredge(max_ml_d, evaluate = TRUE, rank = "AIC")

max_ml_s <- gls(start_lag ~ cum_precip_pre + diurnal_temp_range + 
  evenness + richness + richness:cluster,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")
dredge_s <- dredge(max_ml_s, evaluate = TRUE, rank = "AIC")

max_ml_e <- gls(end_lag ~ cum_precip_end + diurnal_temp_range + 
  evenness + richness + richness:cluster,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")
dredge_e <- dredge(max_ml_e, evaluate = TRUE, rank = "AIC")

max_ml_list <- list(max_ml_c, max_ml_l, max_ml_r, 
  max_ml_d, max_ml_s, max_ml_e)

dredge_list <- list(dredge_c, dredge_l, dredge_r, 
  dredge_d, dredge_s, dredge_e)

# Fit climate only ("null") models
null_ml_c <- gls(cum_vi ~ cum_precip_seas + diurnal_temp_range,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

null_ml_l <- gls(s1_length ~ cum_precip_seas + diurnal_temp_range,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

null_ml_r <- gls(s1_green_rate ~ cum_precip_pre + diurnal_temp_range,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

null_ml_d <- gls(s1_senes_rate ~ cum_precip_end + diurnal_temp_range,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

null_ml_s <- gls(start_lag ~ cum_precip_end + diurnal_temp_range,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

null_ml_e <- gls(end_lag ~ cum_precip_end + diurnal_temp_range,
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

null_ml_list <- list(null_ml_c, null_ml_l, null_ml_r, 
  null_ml_d, null_ml_s, null_ml_e)

# Fit ML version of "best" models
#dredge_list

best_ml_c <- gls(cum_vi ~ cum_precip_seas + diurnal_temp_range + 
  evenness + richness + richness:cluster, 
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

best_ml_l <- gls(s1_length ~ cum_precip_seas + 
  evenness + richness + richness:cluster, 
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

best_ml_r <- gls(s1_green_rate ~ cum_precip_pre + diurnal_temp_range + 
  richness + richness:cluster, 
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

best_ml_d <- gls(s1_senes_rate ~ cum_precip_end + diurnal_temp_range + 
  richness + richness:cluster, 
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

best_ml_s <- gls(start_lag ~ cum_precip_pre + diurnal_temp_range + 
  richness + richness:cluster, 
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

best_ml_e <- gls(end_lag ~ diurnal_temp_range + 
  richness + richness:cluster, 
  correlation = corGaus(1, form = ~x+y),
  data = dat_std, method = "ML")

best_ml_list <- list(best_ml_c, best_ml_l, best_ml_r, 
  best_ml_d, best_ml_s, best_ml_e)

# Fit REML version of "best" models
best_reml_c <- update(best_ml_c, method = "REML")
best_reml_l <- update(best_ml_l, method = "REML")
best_reml_r <- update(best_ml_r, method = "REML")
best_reml_d <- update(best_ml_d, method = "REML")
best_reml_s <- update(best_ml_s, method = "REML")
best_reml_e <- update(best_ml_e, method = "REML")

best_reml_list <- list(best_reml_c, best_reml_l, best_reml_r, 
  best_reml_d, best_reml_s, best_reml_e)

# Get model fit statistics 
modFit <- function(best_ml, null_ml) {
  return(
    data.frame(mod = c("best_ml", "null_ml"),
      aic = c(AIC(best_ml, null_ml)[,2]), 
      bic = c(BIC(best_ml, null_ml)[,2]), 
      rsq = do.call(rbind, lapply(list(best_ml, null_ml), r.squaredLR)),
      logl = c(logLik(best_ml)[1], logLik(null_ml)[1])
    )
  )
}

fit_list <- lapply(seq(length(best_ml_list)), function(x) {
  modFit(best_ml_list[[x]], null_ml_list[[x]])
  })

# Combine all to one list
all_mod_list <- list(max_ml_list, dredge_list, null_ml_list, 
  best_ml_list, best_reml_list, fit_list)

# Write models
saveRDS(all_mod_list, "dat/all_mod_list.rds")
all_mod_list <- readRDS("dat/all_mod_list.rds")

# Export model statistics table
mod_stat_df <- as.data.frame(do.call(rbind, lapply(all_mod_list[[6]], function(x) {
      ##' If positive, max mod better
      daic <- x[2,2] - x[1,2]
      dbic <- x[2,3] - x[1,3]
      rsq <- x[1,4]
      dlogl <- x[2,5] - x[1,5]

  unlist(c(daic, dbic, rsq, dlogl))
}))) 
names(mod_stat_df) <- c("daic", "dbic", "rsq", "dlogl")
mod_stat_df$resp <- unlist(lapply(all_mod_list[[4]], function(x) {
    gsub("\\s~.*", "", x$call[[2]])[2]
  }))
mod_stat_df <- mod_stat_df[,c(length(mod_stat_df), seq(length(mod_stat_df) - 1))]
mod_stat_df$resp <- factor(mod_stat_df$resp, levels = names(resp_lookup), 
  labels = resp_lookup)
mod_stat_tab <- xtable(mod_stat_df, 
  label = "mod_stat",
  align = "rrcccc",
  display = c("s", "s", "f", "f", "f", "f"),
  digits = c(0, 0, 1, 1, 2, 2),
  caption = "Model fit statistics for each phenological metric.")
names(mod_stat_tab) <- c("Response", "$\\delta$AIC", "$\\delta$BIC", "R\\textsuperscript{2}\\textsubscript{adj}", "$\\delta$logLik")

fileConn <- file("out/mod_stat.tex")
writeLines(print(mod_stat_tab, include.rownames = FALSE, 
    table.placement = "H",
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Model selection tables
lapply(all_mod_list[[2]], function(x) {
  out <- as.data.frame(x)[1:10, c(2:5, 8:11)] %>%
    mutate(across(1:4, ~ if_else(is.na(.x), "", "\\checkmark")),
      rank = seq(nrow(.)), .before = 1)

  resp <- gsub("\\s~.*", "", attr(x, "model.calls")[[1]][[2]])[2]

  tab_name <- paste0("mod_sel_", resp)

  out_tab <- xtable(out,
    label = tab_name,
    caption = paste(resp_lookup[names(resp_lookup) %in% resp], 
      "model selection candidate models, with fit statistics."),
    align = "ccccccrrrr",
    display = c("s", "d", "s", "s", "s", "s", "d", "d", "f", "f"),
    digits = c(0,0,0,0,0,0,0,0,2,3))
  
  names(out_tab) <- c("Rank", "Precipitation", "Diurnal dT", "Evenness", 
    "Richness", "logLik", "AIC", "$\\Delta{}IC$", "$W_{i}$")

  fileConn <- file(file.path("out", paste0(tab_name, ".tex")))
  writeLines(print(out_tab, include.rownames = FALSE, 
      sanitize.text.function = function(x) {x}), 
    fileConn)
  close(fileConn)
})

# Combine all tables into one file
writeLines(
  do.call(c, lapply(list.files("out", "^mod_sel.*.tex", full.names = TRUE), readLines)),
  "out/all_mod_sel.tex")

# Model slope plots
mod_slope_df <- do.call(rbind, lapply(all_mod_list[[5]], function(x) {
  mod_pred <- get_model_data(x, type = "est")

  mod_pred_fil <- mod_pred %>%
    filter(!grepl("richness:cluster", term)) %>%
    mutate(resp = gsub("\\s~.*", "", x$call[[2]])[2]) %>%
    dplyr::select(1:8, group, resp)

  return(mod_pred_fil)
  })) %>%
  mutate(term = factor(term, levels = rev(names(pred_lookup)), labels = rev(pred_lookup)),
    resp = factor(resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)]),
    psig = case_when(
      p.value <= 0.05 ~ "*",
      p.value <= 0.01 ~ "**",
      p.value <= 0.001 ~ "***",
      TRUE ~ NA_character_))

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
  labs(x = "Estimate", y = "")
dev.off()

# Export interaction effect on richness plots
rand_marg_df <- do.call(rbind, lapply(all_mod_list[[5]], function(x) {
  preds <- as.data.frame(ggemmeans(model = x,
      terms = c("richness", "cluster")))
  preds$resp <- gsub("\\s~.*", "", x$call[[2]])[2]  
  names(preds) <- c("richness", "pred", "se", "conf_lo", "conf_hi", "cluster", "resp")
  return(preds)
}))
rand_marg_df$resp_clean <- factor(rand_marg_df$resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)])

marg_signif <- do.call(rbind, lapply(all_mod_list[[5]], function(x) {
  summ <- Anova(x)
  psig <- case_when(
      summ[,3] <= 0.05 ~ "*",
      summ[,3] <= 0.01 ~ "**",
      summ[,3] <= 0.001 ~ "***",
      TRUE ~ NA_character_)
  out <- data.frame(resp = gsub("\\s~.*", "", x$call[2]),
    variable = row.names(summ), psig) 
  out <- out[out$variable == "richness:cluster",]
  out$x <- 0
  out$resp
  midpoint <- mean(rand_marg_df[rand_marg_df$resp == out$resp & 
    rand_marg_df$richness == 0,"pred"])
  out$y <- midpoint + 
    (max(rand_marg_df[rand_marg_df$resp == out$resp,"conf_hi"]) - midpoint) * 0.25
  return(out)
  }))
marg_signif$resp_clean <- factor(marg_signif$resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)])

pdf(file = "img/mod_marg.pdf", width = 12, height = 6)
ggplot() + 
  geom_ribbon(data = rand_marg_df, 
    aes(x = richness, ymin = conf_lo, ymax = conf_hi, colour = cluster), 
    linetype = 2, alpha = 0.2, fill = NA) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) +
  new_scale_colour() + 
  geom_line(data = rand_marg_df, 
    aes(x = richness, y = pred, colour = cluster),
    size = 1.5) + 
  geom_text(data = marg_signif, 
    aes(x = x, y = y, label = psig),
    size = 8) + 
  scale_colour_manual(name = "Cluster", 
    values = brightness(clust_pal, 0.75)) +
  facet_wrap(~resp_clean, scales = "free_y") + 
  theme_panel() +
  labs(x = "Species richness", y = "")
dev.off()


