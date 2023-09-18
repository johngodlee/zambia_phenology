# Define predictor name lookup
pred_lookup <- c(
  eff_rich = "Species diversity",
  cluster = "Vegetation type",
  diam_quad_mean = "Mean stem diameter",
  Detarioideae = "Detarioid relative abundance",
  cum_precip_seas = "Rainy season precipitation",
  cum_precip_pre = "Pre-green-up precipitation",
  cum_precip_end = "Pre-senescence precipitation",
  cum_temp_seas = "Rainy season degree days",
  cum_temp_pre = "Pre-green-up degree days",
  cum_temp_end = "Pre-senescence degree days")

# Define response name lookup
resp_lookup <- c(
  EVI_Area = "Cumulative EVI",
  season_length = "Season length",
  green_rate = "Green-up length",
  senes_rate = "Senescence length",
  start_lag = "Pre-rain green-up",
  end_lag = "Senescence lag")

# Define response plot axis names
resp_plot_axes <- resp_lookup
resp_plot_axes[2:6] <- paste(resp_plot_axes[2:6], "(d)")

# Cluster name lookup
clust_lookup <- c(
  "1" = "Uapaca miombo",
  "2" = "Combretaceae woodland",
  "3" = "Julbernardia miombo",
  "4" = "Cryptosepalum miombo")

# Theme colours
pal <- c("lightseagreen", "#DE6400", "dodgerblue", "tomato", "grey", "#E0E0E0")

# Cluster colours
clust_pal <- c(
  "#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#13447b", "#DAA51B", 
  "#2F8AC4", "#764E9F", "#ED645A", "#CC3A8E", "#A5AA99")

# ggplot2 theme
theme_panel <- function() {
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = NA),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14)
  )
}

