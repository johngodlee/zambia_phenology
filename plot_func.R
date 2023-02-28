# Define variable name lookup
pred_lookup <- c(
  eff_rich = "Diversity",
  evenness = "Evenness",
  cum_precip_seas = "Wet season precip.",
  cum_precip_pre = "Pre-green-up precip.",
  cum_precip_end = "Pre-senescence precip.",
  cluster = "Vegetation type",
  diurnal_temp_range = "Diurnal dT",
  diam_quad_mean = "Stem diameter",
  Detarioideae = "Detarioid BA")

resp_lookup <- c(
  EVI_Area = "Cumulative EVI",
  season_length = "Season length",
  green_rate = "Green-up rate",
  senes_rate = "Senescence rate",
  start_lag = "Green-up lag",
  end_lag = "Senescence lag")

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
  )
}


