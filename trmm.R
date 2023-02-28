# Extract TRMM rainy data statistics 
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-02-24

# Packages
library(dplyr)
library(RcppRoll)
library(parallel)
library(english)

source("tex_func.R")

#' Fast rbind for dataframes
#'
#' @param x list of dataframe objects
#'
#' @return dataframe
#' 
fastRbind <- function(x) { 
  list2DF(lapply(setNames( seq_along(x[[1]]), names(x[[1]])), 
      function(i) {
        unlist(lapply(x, `[[`, i), FALSE, FALSE)
      }))
}

# Import data
plots <- readRDS("dat/plots.rds")

# Read raw precipitation time series for each plot 
trmm <- readRDS("dat/trmm_ts.rds")

# Are all plots in the time series data?
stopifnot(all(plots$plot_cluster %in% trmm$plot_cluster))

# Find plots with any NA or negative precipitation values
trmm_any_bad <- trmm %>% 
  filter(is.na(precip)) %>% 
  pull(plot_cluster) %>% 
  unique()

# Check not too many plots are excluded
stopifnot(length(trmm_any_bad) <= 10)

# Exclude plots with dodgy precipitation time series
trmm_fil <- trmm %>% 
  filter(!plot_cluster %in% trmm_any_bad) %>%
  mutate(precip = ifelse(precip < 0, 0, precip))

# Create month and year IDs
trmm_fil$year <- format(trmm_fil$date, "%Y")
trmm_fil$month <- format(trmm_fil$date, "%m")

# Split by plot cluster
trmm_list <- split(trmm_fil, trmm_fil$plot_cluster)

onset_period_one <- 10
onset_precip_one <- 20
onset_period_two <- 20
onset_precip_two <- 20
rainy_def <- 0.5
end_period <- 30
end_rainy_days <- 4

# For each site:
hydro_years_all <- fastRbind(mclapply(trmm_list, function(x) {
  # For each year, find month(s) with least rain
  month_least_rain <- x %>% 
    group_by(year, month) %>%
    summarise(
      precip = sum(precip)) %>% 
    filter(precip == min(precip)) 

  # For each calendar year
  hydro_list <- lapply(seq(2001, 2019, 1), function(y) {
    message(y)
    # Find month with least rain in previous year, define hydro year start date
    prev_year <- month_least_rain[month_least_rain$year == y-1,]
    prev_year_month <- prev_year[prev_year$month == max(prev_year$month),]
    hydro_start <- as.Date(paste(prev_year_month$year, 
        prev_year_month$month, "01", sep = "-"))

    # Find month with least rain in this year, define hydro year end date
    this_year <- month_least_rain[month_least_rain$year == y,]
    this_year_month <- this_year[this_year$month == min(this_year$month),]
    hydro_end <- as.Date(paste(this_year_month$year, 
        as.numeric(this_year_month$month)+1, "01", sep = "-")) - 1
    x[x$date >= hydro_start & x$date <= hydro_end,]
  })

  # For each hydrological year:
  fastRbind(lapply(hydro_list, function(y) { 
    # Calculate 10 day rolling sum of precipitation
    y$roll10 <- roll_suml(y$precip, onset_period_one)

    # Calculate rolling sum of rainy days
    y$rainy30 <- roll_suml(y$precip >= rainy_def, end_period)

    # Find dates with >10 days with >25 mm rainfall
    day10_25 <- y[y$roll10 > onset_precip_one, "date"]

    y$day10_25 <- y$date %in% day10_25

    # Find first subsequent 20 day period with >20 mm rainfall
    good_start <- FALSE
    rep_date <- 1
    while (good_start == FALSE) {
      search_start <- day10_25[rep_date] + 10
      search_end <- search_start + onset_period_two
      precip_sum <- sum(y[y$date > search_start & y$date <= search_end, "precip"])
      if (precip_sum > onset_precip_two) {
        good_start <- TRUE
        trmm_start <- day10_25[rep_date]
      } else {
        rep_date <- rep_date + 1
      }
    }

    # Filter to after start of rainy season (+90 day buffer)
    y_fil <- y[y$date > trmm_start + 90,]

    # Find first date after start of rainy season with <4 rainy days in a 30 day period.
    trmm_end <- min(y_fil[y_fil$rainy30 < end_rainy_days, "date"], 
      na.rm = TRUE)

    # Filter to hydrological year
    y_hydro <- y_fil[y_fil$date < trmm_end,]

    # Calculate cumulative rainfall
    y_hydro$cum_rain <- cumsum(y_hydro$precip)

    # Find dates where cumulative rainfall >10% and >95% of total
    trmm_start10 <- y_hydro[y_hydro$cum_rain > 0.1*sum(y_hydro$precip), "date"][1]
    trmm_end95 <- y_hydro[y_hydro$cum_rain > 0.95*sum(y_hydro$precip), "date"][1]

    # Calculate precipitation in 90 days before start and end of rainy season
    cum_precip_pre <- sum(y[y$date > trmm_start-90 & y$date < trmm_start, "precip"])
    cum_precip_end <- sum(y[y$date > trmm_end-90 & y$date < trmm_end, "precip"])

    # Create pretty dataframe
    data.frame(
      plot_cluster = unique(y$plot_cluster),
      year = max(y$year),
      trmm_start,
      trmm_end,
      trmm_start10,
      trmm_end95,
      trmm_length = trmm_end - trmm_start,
      cum_precip_seas = sum(y_hydro$precip),
      cum_precip_pre,
      cum_precip_end)
    
  }))
}, mc.cores = detectCores()-1))

# Filter to good plots and fix dates
hydro_years_all_fil <- hydro_years_all %>% 
  filter(plot_cluster %in% plots$plot_cluster) %>% 
  mutate(
    across(contains(c("trmm_start", "trmm_end")), 
      ~as.numeric(as.Date(.x, origin = "1970-01-01") - 
        as.Date(paste0(year, "-01-01")))))

# LaTeX variables
write(
  c(
    commandOutput(onset_period_one, "onsetPeriodOne"),
    commandOutput(onset_precip_one, "onsetPrecipOne"),
    commandOutput(onset_period_two, "onsetPeriodTwo"),
    commandOutput(onset_precip_two, "onsetPrecipTwo"),
    commandOutput(rainy_def, "rainyDef"),
    commandOutput(end_period, "endPeriod"),
    commandOutput(as.english(end_rainy_days), "endRainyDays")
  ),
  file = "out/trmm_vars.tex")

# Write data 
saveRDS(hydro_years_all_fil, "dat/trmm.rds")

