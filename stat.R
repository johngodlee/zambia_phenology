# Extract rainy season and growing season statistics
# John Godlee (johngodlee@gmail.com)
# Last updated: 2023-02-24

# Packages
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(RcppRoll)
library(parallel)
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
era <- readRDS("./dat/era_ts.rds")
trmm <- readRDS("dat/trmm_ts.rds")
modis <- readRDS("./dat/modis_ts.rds")

# Create a value for Dormancy which is just MidGreendown + 30 days
dorm_pad <- 30
##' Because mean diff between greenup and mid-greenup is 35.42+/-6.379 days
modis$Dormancy2 <- modis$MidGreendown + dorm_pad

# Calculate length of greening and senescence period
modis$green_rate <- modis$Maturity - modis$Greenup
modis$senes_rate <- modis$Dormancy2 - modis$Senescence
modis$season_length <- modis$Dormancy2 - modis$Greenup

# Check all scenese of "best" quality
stopifnot(all(modis$QA_Overall == 0 | is.na(modis$QA_Overall)))

# Get range of data
modis_range <- range(modis$year)

# Are all plots in the time series data?
stopifnot(all(plots$plot_cluster %in% trmm$plot_cluster))

# Find plots with any NA values
trmm_any_bad <- trmm %>% 
  filter(is.na(precip)) %>% 
  pull(plot_cluster) %>% 
  unique()

# Check not too many plots are excluded
stopifnot(length(trmm_any_bad) <= 10)

# Exclude plots with dodgy precipitation time series
trmm_fil <- trmm %>% 
  filter(
    !plot_cluster %in% trmm_any_bad,
    plot_cluster %in% plots$plot_cluster) %>% 
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

# Turn off messages from dplyr::summarise()
options(dplyr.summarise.inform = FALSE)

# For each site:
stat_all <- fastRbind(mclapply(trmm_list, function(x) {
  # Display site ID
  message(unique(x$plot_cluster))

  # For each year, find month(s) with least rain
  month_least_rain <- x %>% 
    group_by(year, month) %>%
    summarise(precip = sum(precip)) %>% 
    filter(precip == min(precip)) 

  # For each calendar year:
  hydro_list <- lapply(seq(2001, 2019, 1), function(y) {
    # Find month with least rain in previous year
    prev_year <- month_least_rain[month_least_rain$year == y-1,]
    prev_year_month <- prev_year[prev_year$month == max(prev_year$month),]

    # define hydrological year start date
    hydro_start <- as.Date(paste(prev_year_month$year, 
        prev_year_month$month, "01", sep = "-"))

    # Find month after with least rain
    this_year <- month_least_rain[month_least_rain$year == y,]
    this_year_month <- this_year[this_year$month == min(this_year$month),]

    # define hydrological year end date
    hydro_end <- as.Date(paste(this_year_month$year, 
        as.numeric(this_year_month$month)+1, "01", sep = "-")) - 1

    # Return TRMM data in hydrological years
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

    # Define rainy season end date
    # First date after start of rainy season (+90 day buffer),
    # with <4 rainy days in a 30 day period.
    trmm_end <- min(y[y$rainy30 < end_rainy_days & 
      y$date > trmm_start + 90, "date"], na.rm = TRUE)

    # Filter to rainy season
    y_rainy <- y[y$date <= trmm_end & y$date > trmm_start,]

    # Calculate cumulative rainfall in hydrological year
    y$cum_rain <- cumsum(y$precip)

    # Find dates where cumulative rainfall >10% and >95% of total
    trmm_start10 <- y[y$cum_rain > 0.1 * sum(y$precip), "date"][1]
    trmm_end95 <- y[y$cum_rain > 0.95 * sum(y$precip), "date"][1]

    # Get MODIS data
    modis_fil <- modis[
      modis$plot_cluster == unique(y_rainy$plot_cluster) & 
        modis$year == max(y_rainy$year),]

    # Format greenup and senescence dates as Dates
    greenup_date <- as.Date(modis_fil$Greenup, origin = "1970-01-01")
    senes_date <- as.Date(modis_fil$Senescence, origin = "1970-01-01")

    # Cumulative precipitation in 30 days before start and end of growing season
    cum_precip_pre <- sum(y[
      y$date >= greenup_date-30 & y$date < greenup_date, "precip"])
    cum_precip_end <- sum(y[
      y$date >= senes_date-30 & y$date < senes_date, "precip"])

    # Subset temperature data
    era_fil <- era[era$plot_cluster == unique(y_rainy$plot_cluster),]

    # Cumulative temperature in 30 days before start and end of growing season 
    cum_temp_pre <- sum(era_fil[
      era_fil$date >= greenup_date-30 & 
      era_fil$date < greenup_date, "era_temp"], na.rm = TRUE)
    cum_temp_end <- sum(era_fil[
      era_fil$date >= senes_date-30 & 
      era_fil$date < senes_date, "era_temp"], na.rm = TRUE)
    cum_temp_seas <- sum(era_fil[
      era_fil$date >= greenup_date & 
      era_fil$date < senes_date, "era_temp"], na.rm = TRUE)

    # Create pretty dataframe
    out <- data.frame(
      plot_cluster = unique(y$plot_cluster),
      year = max(y$year),
      trmm_start,
      trmm_end,
      trmm_start10,
      trmm_end95,
      trmm_length = trmm_end - trmm_start,
      cum_precip_seas = sum(y_rainy$precip),
      cum_precip_pre,
      cum_precip_end,
      cum_temp_seas,
      cum_temp_pre,
      cum_temp_end)

    cbind(out, modis_fil[,!names(modis_fil) %in% c("plot_cluster", "year")])
  }))
}, mc.cores = detectCores()-1))

# Define date variables
date_vars <- c("Greenup", "MidGreenup", "Maturity", "Senescence", 
  "MidGreendown", "Dormancy2")

# Fix dates
# Calculate lags
stat_all_clean <- stat_all %>% 
  mutate(
    across(contains(c("trmm_start", "trmm_end")), 
      ~as.numeric(as.Date(.x, origin = "1970-01-01") - 
        as.Date(paste0(year, "-01-01")))),
    cum_temp_pre = cum_temp_pre / 30,
    cum_temp_end = cum_temp_end / 30,
    cum_temp_seas = cum_temp_seas / season_length,
    across(all_of(date_vars), 
      ~as.Date(.x, origin = "1970-01-01") - as.Date(paste0(year, "-01-01"))),
    across(all_of(date_vars), 
      ~ifelse(.x > 182, .x-365, .x)),
    start_lag = -(Greenup - trmm_start),
    end_lag = Dormancy2 - trmm_end) %>%
  filter(end_lag > -100)

# LaTeX variables
write(
  c(
    commandOutput(onset_period_one, "onsetPeriodOne"),
    commandOutput(onset_precip_one, "onsetPrecipOne"),
    commandOutput(onset_period_two, "onsetPeriodTwo"),
    commandOutput(onset_precip_two, "onsetPrecipTwo"),
    commandOutput(rainy_def, "rainyDef"),
    commandOutput(end_period, "periodEnd"),
    commandOutput(end_rainy_days, "rainyDaysEnd"),
    commandOutput(modis_range[1], "modisStart"),
    commandOutput(modis_range[2], "modisEnd")
  ),
  file = "out/stat_vars.tex")

# Write data 
saveRDS(stat_all_clean, "dat/stat_all.rds")
##' For each year for each plot,
##' start and end of rainy season, defined using rainy days
##' start and end of rainy season, defined using quantiles of annual precip.
##' rainy season length
##' Cumulative precip. before season, during season, end of season
##' Cumulative temp. (degree days) before season, during season, end of season

# Calculate average values and lag times
out_avg <- stat_all_clean %>% 
  dplyr::select(-QA_Overall, -year) %>%
  group_by(plot_cluster) %>% 
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

# Write extracted values to file
saveRDS(out_avg, "dat/stat_avg.rds")

