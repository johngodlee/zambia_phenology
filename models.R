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
precip_vars <- c(
  "cum_precip_seas", 
  "cum_precip_seas", 
  "cum_precip_pre",
  "cum_precip_end", 
  "cum_precip_pre", 
  "cum_precip_seas")

other_vars <- "diurnal_temp_range + evenness + eff_rich + diam_quad_mean + evenness:cluster + eff_rich:cluster + diam_quad_mean:cluster + cluster"

max_mod_flist <- paste0(names(resp_lookup), " ~ ", precip_vars, " + ", other_vars)

max_ml_list <- lapply(max_mod_flist, function(x) {
  lm(x, data = dat_std, na.action = na.fail)
  })

dredge_list <- lapply(max_ml_list, function(x) {
  dredge(x, evaluate = TRUE, rank = "AIC")
  })

# Fit climate only "null" models
null_mod_flist <- paste0(names(resp_lookup), " ~ ", precip_vars, " + ", "diurnal_temp_range")

null_ml_list <- lapply(null_mod_flist, function(x) {
  lm(x, data = dat_std)
  })

# Fit REML version of "best" models
##' For marginal effect slopes
# Which models are "best" while including interaction of cluster and richness?

best_mod_flist <- c(
  "cum_vi ~ cum_precip_seas + evenness + eff_rich + diam_quad_mean + evenness:cluster + diam_quad_mean:cluster + cluster",
  "s1_length ~ cum_precip_seas + evenness + eff_rich + diam_quad_mean + evenness:cluster + cluster",
  "s1_green_rate ~ diurnal_temp_range + cluster",
  "s1_senes_rate ~ cum_precip_end + cluster",
  "start_lag ~ cum_precip_pre + diurnal_temp_range + eff_rich + evenness + cluster",
  "end_lag ~ cum_precip_seas + diurnal_temp_range + diam_quad_mean + eff_rich + diam_quad_mean + eff_rich:cluster + cluster")

best_ml_list <- lapply(best_mod_flist, function(x) {
  lm(x, data = dat_std)
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
  modFit(max_ml_list[[x]], null_ml_list[[x]])
  })

# Nest all model lists in one list
all_mod_list <- list(max_ml_list, dredge_list, null_ml_list, 
  best_ml_list, fit_list)
names(all_mod_list) <- c("max_ml", "dredge", "null_ml",
  "best_ml", "fit_list")

# Write model lists
saveRDS(all_mod_list, "dat/all_mod_list.rds")

# Export model statistics table
mod_stat_df <- as.data.frame(do.call(rbind, lapply(all_mod_list[[5]], function(x) {
      # If positive, max mod better
      daic <- x[2,2] - x[1,2]
      dbic <- x[2,3] - x[1,3]
      rsq <- x[1,4]
      dlogl <- x[2,5] - x[1,5]

  unlist(c(daic, dbic, rsq, dlogl))
}))) 

names(mod_stat_df) <- c("daic", "dbic", "rsq", "dlogl")
mod_stat_df$resp <- names(resp_lookup)
mod_stat_df <- mod_stat_df[,c(5,1,3,4)]
mod_stat_df$resp <- factor(mod_stat_df$resp, levels = names(resp_lookup), 
  labels = resp_lookup)
mod_stat_tab <- xtable(mod_stat_df, 
  label = "mod_stat",
  align = "rrccc",
  display = c("s", "s", "f", "f", "f"),
  digits = c(0, 0, 1, 2, 2),
  caption = "Model fit statistics for the best model describing each phenological metric.")
names(mod_stat_tab) <- c("Response", "$\\delta$AIC", "R\\textsuperscript{2}\\textsubscript{adj}", "$\\delta$logLik")

fileConn <- file("out/mod_stat.tex")
writeLines(print(mod_stat_tab, include.rownames = FALSE, 
    table.placement = "H",
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Model selection tables
# Highlight best model according to AIC
best_mod <- c(1,2,1,2,1,2)  # best_ml_list

lapply(seq(length(all_mod_list[[2]])), function(x) {
  out <- as.data.frame(all_mod_list[[2]][[x]])[1:10, c(2:15)] %>%
    mutate(across(1:9, ~ if_else(is.na(.x), "", "\\checkmark")),
      rank = seq(nrow(.)), .before = 1) %>%
    mutate(
      diam_quad_mean = if_else(`cluster:diam_quad_mean` == "\\checkmark", "\\checkmark+", diam_quad_mean),
      eff_rich = if_else(`cluster:eff_rich` == "\\checkmark", "\\checkmark+", eff_rich),
      evenness = if_else(`cluster:evenness` == "\\checkmark", "\\checkmark+", evenness)) %>%
    dplyr::select(-starts_with("cluster:"))

  names(out)[3] <- "precip"

  out <- out %>%
    mutate(rank = sprintf("%.0f", rank),
      logLik = sprintf("%.0f", logLik),
      AIC = sprintf("%.0f", AIC),
      delta = sprintf("%.0f", delta), 
      weight = sprintf("%.3f", weight))

  resp <- gsub("\\s~.*", "", attr(all_mod_list[[2]][[x]], "model.calls")[[1]][[2]])[2]

  tab_name <- paste0("mod_sel_", resp)

  out[best_mod[x],] <- paste0("\\underline{\\textbf{", out[best_mod[x],], "}}")

  out_tab <- xtable(out,
    label = tab_name,
    caption = paste(resp_lookup[names(resp_lookup) %in% resp], 
      "model selection candidate models, with fit statistics. The overall best model is marked by bold text, according to AIC and model parsimony."),
    align = "cccccccccrrrr",
    display = c("s", "d", "s", "s", "s", "s", "s", "s", "s", "d", "d", "f", "f"),
    digits = c(0,0,0,0,0,0,0,0,0,0,0,2,3))
  
  names(out_tab) <- c("Rank", "Cluster", "Precipitation", "Stem diameter", "Diurnal dT", "Richness", "Evenness", 
    "DoF", "logLik", "AIC", "$\\Delta{}$", "$W_{i}$")

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
intf <- function(x) {
  best_int_mod_flist <- best_mod_flist[grepl(x, best_mod_flist)] %>%
    ifelse(grepl(paste0(x, ":cluster"), .), ., paste0(., " + ", x, ":cluster"))

  best_int_ml_list <- lapply(best_int_mod_flist, function(x) {
    lm(formula(x), data = dat_std)
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
evenness_intf <- intf("evenness")
diam_quad_mean_intf <- intf("diam_quad_mean")

marg_df <- rbind(eff_rich_intf[[2]], evenness_intf[[2]], diam_quad_mean_intf[[2]])

marg_df$resp_clean <- factor(marg_df$resp, 
      levels = names(resp_lookup[c(1,3,5,2,4,6)]), 
      labels = resp_lookup[c(1,3,5,2,4,6)])

marg_df$pred_name_clean <- factor(marg_df$pred_name, 
      levels = names(pred_lookup), 
      labels = pred_lookup)

pdf(file = "img/mod_marg.pdf", width = 8, height = 10)
ggplot() + 
  geom_ribbon(data = marg_df, 
    aes(x = val, ymin = conf_lo, ymax = conf_hi, colour = cluster), 
    linetype = 2, alpha = 0.2, fill = NA) + 
  scale_colour_manual(name = "Cluster", values = clust_pal) +
  new_scale_colour() + 
  geom_line(data = marg_df, 
    aes(x = val, y = pred, colour = cluster),
    size = 1.5) + 
  scale_colour_manual(name = "Cluster", 
    values = brightness(clust_pal, 0.75)) +
  facet_grid(resp_clean~pred_name_clean, scales = "free") + 
  theme_panel() + 
  labs(x = "", y = "")
dev.off()

# Post-hoc tests for significance of slope differences among clusters
rich_lsq_list <- lapply(eff_rich_intf[[1]], function(x) {
  emmeans(x, 
  pairwise ~ cluster*eff_rich, 
  adjust = "tukey")
})

# Extract terms
lsq_terms_df <- do.call(rbind, lapply(rich_lsq_list, function(x) {
  out <- as.data.frame(x[[2]])
  out$resp <- names(attributes(x$emmeans@model.info$terms)$dataClasses[1])
  return(out)
}))

# Make tidy
lsq_terms <- lsq_terms_df %>%
  dplyr::select(resp, contrast, estimate, SE, df, t.ratio, p.value) %>%
  mutate(p.value = pFormat(p.value))

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

# Create table
lsq_terms_tab <- xtable(lsq_terms, 
  label = "lsq_terms",
  align = "rrcccccc",
  display = c("s", "s", "s", "E", "E", "d", "f", "f"),
  digits = c(  0,   0,   0,   1,   2,   0,   2,   2 ),
  caption = "Comparisons of species diversity interaction marginal effects using post-hoc Tukey's tests.")
names(lsq_terms_tab) <- c("Response", "Clusters", "Estimate", "SE", "DoF", 
  "T ratio", "Prob.")

# Write table to file
fileConn <- file("out/lsq_terms.tex")
writeLines(print(lsq_terms_tab, include.rownames = FALSE, 
    table.placement = "H",
    hline.after = c(-1,0,
      seq(from = clust_combn, 
        by = clust_combn, 
        length.out = length(unique(lsq_terms$resp)) - 1)),
    sanitize.text.function = function(x) {x}), 
  fileConn)
close(fileConn)

# Write variables
write(
  c(
    commandOutput(floor(mod_stat_df$rsq[mod_stat_df$resp == "Green-up lag"] * 100), "greenLagrsq"),
    commandOutput(floor(mod_stat_df$rsq[mod_stat_df$resp == "Senescence lag"] * 100), "senesLagrsq")
    ),
  file = "out/models_vars.tex")

