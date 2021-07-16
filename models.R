# Testing effect of species composition and richness on land surface phenology with Zambia ILUAii data 
# John Godlee (johngodlee@gmail.com)
# 2020-10-25

# Set working directory

# Packages
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(ggeffects)
library(shades)
library(xtable)
library(lme4)
library(MuMIn)
library(sjPlot)
library(car)
library(emmeans)

source("functions.R")

# Import data
dat_std <- readRDS("dat/plots_anal.rds")

# Mixed models of diversity with clusters 

# Fit maximal models and rank model subsets by AIC
max_mod_flist <- paste0(names(resp_lookup), " ~ cum_precip_seas + diurnal_temp_range + evenness + eff_rich + (1|cluster)")

max_ml_list <- lapply(max_mod_flist, function(x) {
  lmer(x, data = dat_std, REML = FALSE, na.action = na.fail)
  })

dredge_list <- lapply(max_ml_list, function(x) {
  dredge(x, evaluate = TRUE, rank = "AIC")
  })

# Fit climate only "null" models
null_mod_flist <- paste0(names(resp_lookup), " ~ cum_precip_seas + diurnal_temp_range + (1|cluster)")

null_ml_list <- lapply(null_mod_flist, function(x) {
  lmer(x, data = dat_std, REML = FALSE)
  })

# Fit REML version of "best" models
##' For marginal effect slopes
# Which models are "best" while including interaction of cluster and richness?

best_mod_flist <- c(
  "cum_vi ~ cum_precip_seas + evenness + eff_rich + (1|cluster)",
  "s1_length ~ cum_precip_seas + diurnal_temp_range + evenness + eff_rich + (1|cluster)",
  "s1_green_rate ~ cum_precip_pre + diurnal_temp_range + (1|cluster)",
  "s1_senes_rate ~ cum_precip_end + evenness + eff_rich + (1|cluster)",
  "start_lag ~ cum_precip_pre + diurnal_temp_range + evenness + eff_rich + (1|cluster)",
  "end_lag ~ cum_precip_seas + diurnal_temp_range + (1|cluster)")

best_ml_list <- lapply(best_mod_flist, function(x) {
  lmer(formula(x), data = dat_std, REML = FALSE)
  })

# Fit REML version of "best" models
##' For model slope interval plots
best_reml_list <- lapply(best_ml_list, update, REML = TRUE)

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
  modFit(max_ml_list[[x]], null_ml_list[[x]])
  })

# Nest all model lists in one list
all_mod_list <- list(max_ml_list, dredge_list, null_ml_list, 
  best_ml_list, best_reml_list, fit_list)
names(all_mod_list) <- c("max_ml", "dredge", "null_ml",
  "best_ml", "best_reml", "fit_list")

# Write model lists
saveRDS(all_mod_list, "dat/all_mod_list.rds")

# Export model statistics table
mod_stat_df <- as.data.frame(do.call(rbind, lapply(all_mod_list[[6]], function(x) {
      # If positive, max mod better
      daic <- x[2,2] - x[1,2]
      dbic <- x[2,3] - x[1,3]
      rsq <- x[1,4]
      dlogl <- x[2,5] - x[1,5]

  unlist(c(daic, dbic, rsq, dlogl))
}))) 

names(mod_stat_df) <- c("daic", "dbic", "rsq", "dlogl")
mod_stat_df$resp <- unlist(lapply(all_mod_list[[4]], function(x) {
    gsub("\\s~.*", "", x@call[[2]])[2]
  }))
mod_stat_df <- mod_stat_df[,c(length(mod_stat_df), seq(length(mod_stat_df) - 1))]
mod_stat_df$resp <- factor(mod_stat_df$resp, levels = names(resp_lookup), 
  labels = resp_lookup)
mod_stat_tab <- xtable(mod_stat_df, 
  label = "mod_stat",
  align = "rrcccc",
  display = c("s", "s", "f", "f", "f", "f"),
  digits = c(0, 0, 1, 1, 2, 2),
  caption = "Model fit statistics for the best model describing each phenological metric.")
names(mod_stat_tab) <- c("Response", "$\\delta$AIC", "$\\delta$BIC", "R\\textsuperscript{2}\\textsubscript{adj}", "$\\delta$logLik")

fileConn <- file("out/mod_stat.tex")
writeLines(print(mod_stat_tab, include.rownames = FALSE, 
    table.placement = "H",
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Model selection tables
# Highlight best model according to AIC
best_mod <- c(1,2,1,2,1,1)  # best_ml_list

lapply(seq(length(all_mod_list[[2]])), function(x) {
  out <- as.data.frame(all_mod_list[[2]][[x]])[1:10, c(2:10)] %>%
    mutate(across(1:4, ~ if_else(is.na(.x), "", "\\checkmark")),
      rank = seq(nrow(.)), .before = 1)

  names(out)[2] <- "precip"

  out <- out %>%
    mutate(rank = sprintf("%.0f", rank),
      logLik = sprintf("%.0f", logLik),
      AIC = sprintf("%.0f", AIC),
      delta = sprintf("%.0f", delta), 
      weight = sprintf("%.3f", weight))

  resp <- gsub("\\s~.*", "", attr(all_mod_list[[2]][[x]], "model.calls")[[1]][[2]])[2]

  tab_name <- paste0("mod_sel_", resp)

  out[best_mod[x],] <- paste0("\\textbf{", out[best_mod[x],], "}")

  out_tab <- xtable(out,
    label = tab_name,
    caption = paste(resp_lookup[names(resp_lookup) %in% resp], 
      "Model selection candidate models, with fit statistics. The overall best model is marked by bold text, according to AIC and model parsimony"),
    align = "cccccccrrrr",
    display = c("s", "d", "s", "s", "s", "s", "s", "d", "d", "f", "f"),
    digits = c(0,0,0,0,0,0,0,0,0,2,3))
  
  names(out_tab) <- c("Rank", "Precipitation", "Diurnal dT", "Evenness", 
    "Richness", "Richness:Cluster", "logLik", "AIC", "$\\Delta{}IC$", "$W_{i}$")

  fileConn <- file(file.path("out", paste0(tab_name, ".tex")))
  writeLines(print(out_tab, include.rownames = FALSE, 
    table.placement = "H",
    sanitize.text.function = function(x) {x}), 
    fileConn)
  close(fileConn)
})

# Combine all tables into one file
mod_sel_files <- list.files("out", "^mod_sel.*.tex", full.names = TRUE)
mod_sel_files <- mod_sel_files[c(1,4,3,5,6,2)]
writeLines(
  do.call(c, lapply(mod_sel_files, readLines)),
  "out/all_mod_sel.tex")

# Model slope plots
mod_slope_df <- do.call(rbind, lapply(all_mod_list[[5]], function(x) {
  mod_pred <- get_model_data(x, type = "est")

  mod_pred_fil <- mod_pred %>%
    filter(!grepl("eff_rich:cluster", term)) %>%
    mutate(resp = gsub("\\s~.*", "", x@call[[2]])[2]) %>%
    dplyr::select(1:8, group, resp)

  return(mod_pred_fil)
  })) %>%
  mutate(term = factor(term, levels = rev(names(pred_lookup)), labels = rev(pred_lookup)),
    resp = factor(resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)]),
    psig = case_when(
      conf.low >= 0 & conf.high >= 0 | conf.low <= 0 & conf.high <= 0 ~ "*",
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
best_int_mod_flist <- c(
  "cum_vi ~ cum_precip_seas + evenness + eff_rich + eff_rich:cluster + (1|cluster)",
  "s1_length ~ cum_precip_seas + diurnal_temp_range + evenness + eff_rich + eff_rich:cluster + (1|cluster)",
  "s1_green_rate ~ cum_precip_pre + diurnal_temp_range + eff_rich + eff_rich:cluster + (1|cluster)",
  "s1_senes_rate ~ cum_precip_end + evenness + eff_rich + eff_rich:cluster + (1|cluster)",
  "start_lag ~ cum_precip_pre + diurnal_temp_range + evenness + eff_rich + eff_rich:cluster + (1|cluster)",
  "end_lag ~ cum_precip_seas + diurnal_temp_range + eff_rich + eff_rich:cluster + (1|cluster)")

best_int_ml_list <- lapply(best_int_mod_flist, function(x) {
  lmer(formula(x), data = dat_std, REML = FALSE)
  })

rand_marg_df <- do.call(rbind, lapply(best_int_ml_list, function(x) {
  preds <- ggemmeans(x, terms = c("eff_rich", "cluster"), 
    type = "re", interval = "confidence", ci.lvl = 0.5)
  preds$resp <- gsub("\\s~.*", "", x@call[[2]])[2]  
  names(preds) <- c("eff_rich", "pred", "se", "conf_lo", "conf_hi", "cluster", "resp")
  return(preds)
}))

rand_marg_df$resp_clean <- factor(rand_marg_df$resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)])

pdf(file = "img/mod_marg.pdf", width = 12, height = 6)
ggplot() + 
  geom_ribbon(data = rand_marg_df, 
    aes(x = eff_rich, ymin = conf_lo, ymax = conf_hi, colour = cluster), 
    linetype = 2, alpha = 0.2, fill = NA) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) +
  new_scale_colour() + 
  geom_line(data = rand_marg_df, 
    aes(x = eff_rich, y = pred, colour = cluster),
    size = 1.5) + 
  scale_colour_manual(name = "Cluster", 
    values = brightness(clust_pal, 0.75)) +
  facet_wrap(~resp_clean, scales = "free_y") + 
  theme_panel() +
  labs(x = "Species richness", y = "")
dev.off()

# Post-hoc tests for significance of slope differences among clusters
lsq_list <- lapply(best_int_ml_list, function(x) {
  emmeans(x, 
  pairwise ~ cluster*eff_rich, 
  adjust = "tukey")
})

# Extract terms
lsq_terms_df <- do.call(rbind, lapply(lsq_list, function(x) {
  out <- as.data.frame(x[[2]])
  out$resp <- gsub("\\~.*", "", x[[1]]@model.info$call[[2]])[2]
  return(out)
}))

# Make tidy
lsq_terms <- lsq_terms_df %>%
  dplyr::select(resp, contrast, estimate, SE, df, t.ratio, p.value)

lsq_terms$contrast <- unlist(lapply(
    lapply(strsplit(lsq_terms$contrast, split = " "), `[`, c(1,4)), 
    function(x) {
  gsub("\\(", "", paste(x[1], x[2], sep = "-"))
}))

lsq_terms$resp <- as.character(factor(lsq_terms$resp, levels = names(resp_lookup), labels = resp_lookup))

clust_combn <- dim(combn(seq_along(unique(dat_std$cluster)), 2))[2]

resp_blanks <- seq_len(nrow(lsq_terms))[-which(seq_len(nrow(lsq_terms)) %in% 
  seq(from = floor(median(seq_along(unique(lsq_terms$resp)))), 
  by = clust_combn, length.out = clust_combn))]

lsq_terms[resp_blanks, "resp"] <- ""

lsq_terms_tab <- xtable(lsq_terms, 
  label = "lsq_terms",
  align = "rrcccccc",
  display = c("s", "s", "s", "E", "E", "d", "f", "f"),
  digits = c(  0,   0,   0,   1,   2,   0,   2,   2 ),
  caption = "Comparisons of interaction marginal effects using post-hoc Tukey's tests.")
names(lsq_terms_tab) <- c("Response", "Clusters", "Estimate", "SE", "DoF", 
  "T ratio", "Prob.")

fileConn <- file("out/lsq_terms.tex")
writeLines(print(lsq_terms_tab, include.rownames = FALSE, 
    table.placement = "H",
    hline.after = c(-1,0,
      seq(from = clust_combn, 
        by = clust_combn, 
        length.out = length(unique(resp_lookup)))),
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Models for demographic structure

diam_mod_flist <- paste0(names(resp_lookup), " ~ diam_quad_mean + (1|cluster)")

diam_ml_list <- lapply(diam_mod_flist, function(x) {
  lmer(x, data = dat_std, REML = FALSE)
  })
 
diam_mod_slope_df <- do.call(rbind, lapply(diam_ml_list, function(x) {
  mod_pred <- get_model_data(x, type = "std")

  mod_pred_fil <- mod_pred %>%
    mutate(resp = gsub("\\s~.*", "", x@call[[2]])[2]) %>%
    dplyr::select(1:8, group, resp)

  return(mod_pred_fil)
  })) %>%
  mutate(term = factor(term, levels = rev(names(pred_lookup)), labels = rev(pred_lookup)),
    resp = factor(resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)]),
    psig = case_when(
      conf.low >= 0 & conf.high >= 0 | conf.low <= 0 & conf.high <= 0 ~ "*",
      TRUE ~ NA_character_))

pdf(file = "img/diam_quad_mod_slopes.pdf", width = 6, height = 5)
ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(data = diam_mod_slope_df, 
    aes(xmin = conf.low, xmax = conf.high, y = resp, colour = group),
    height = 0) + 
  geom_point(data = diam_mod_slope_df,
    aes(x = estimate, y = resp, fill = group),
    shape = 21, colour = "black") + 
  geom_text(data = diam_mod_slope_df,
    aes(x = estimate, y = resp, colour = group, label = psig),
    size = 8, nudge_y = 0.1) + 
  scale_fill_manual(name = "", values = pal[3:4]) +
  scale_colour_manual(name = "", values = pal[3:4]) +
  theme_panel() + 
  theme(legend.position = "none") + 
  labs(x = "Estimate", y = "")
dev.off()

diam_int_mod_flist <- paste0(names(resp_lookup), " ~ diam_quad_mean + (diam_quad_mean|cluster)")

diam_int_ml_list <- lapply(diam_int_mod_flist, function(x) {
  lmer(x, data = dat_std, REML = FALSE)
  })

rand_marg_df <- do.call(rbind, lapply(diam_int_ml_list, function(x) {
  preds <- ggpredict(x, terms = c("diam_quad_mean", "cluster"), 
    type = "re", interval = "confidence", ci.lvl = 0.5)
  preds$resp <- gsub("\\s~.*", "", x@call[[2]])[2]  
  names(preds) <- c("diam_quad_mean", "pred", "se", "conf_lo", "conf_hi", "cluster", "resp")
  return(preds)
}))

rand_marg_df$resp_clean <- factor(rand_marg_df$resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)])

pdf(file = "img/diam_mod_marg.pdf", width = 12, height = 6)
ggplot() + 
  geom_ribbon(data = rand_marg_df, 
    aes(x = diam_quad_mean, ymin = conf_lo, ymax = conf_hi, colour = cluster), 
    linetype = 2, alpha = 0.2, fill = NA) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) +
  new_scale_colour() + 
  geom_line(data = rand_marg_df, 
    aes(x = diam_quad_mean, y = pred, colour = cluster),
    size = 1.5) + 
  scale_colour_manual(name = "Cluster", 
    values = brightness(clust_pal, 0.75)) +
  facet_wrap(~resp_clean, scales = "free_y") + 
  theme_panel() +
  labs(x = "Tree size", y = "")
dev.off()
