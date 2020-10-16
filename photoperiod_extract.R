# Extract photoperiod from dataset for each plot
# John Godlee (johngodlee@gmail.com)
# 2020-10-16

# Set working directory

# Packages
library(dplyr)


# Import data
dat <- readRDS("dat/plots_div.rds")

pp <- readRDS("dat/photoperiod_extract.rds")


# Centre s1_start on previous year
dat$s1_start_centre <- 365 - abs(dat$s1_start)

# For each plot, extract cumulative photoperiod 30 days before start of growing season
pp_split <- split(pp, pp$plot_cluster)

prior <- lapply(dat$plot_cluster, function(x) {
  s1_start <- dat[dat$plot_cluster == x, "s1_start_centre"]
  pp_fil <- pp_split[names(pp_split) == x][[1]]
  pp_fil[pp_fil$doy < s1_start &
    pp_fil$doy > s1_start - 30,]
})

pp_prior <- unlist(lapply(prior, function(x) {
  sum(x$photoperiod, na.rm = TRUE)
}))

# For each plot, extract cumulative photoperiod 30 days before start of growing season
after <- lapply(dat$plot_cluster, function(x) {
  s1_end <- dat[dat$plot_cluster == x, "s1_end"]
  pp_fil <- pp_split[names(pp_split) == x][[1]]
  pp_fil[pp_fil$doy < s1_end &
    pp_fil$doy > s1_end - 30,]
})

pp_prior <- unlist(lapply(prior, function(x) {
  sum(x$photoperiod, na.rm = TRUE)
}))
pp_after <- unlist(lapply(after, function(x) {
  sum(x$photoperiod, na.rm = TRUE)
}))

# Create clean dataframe for output
dat_clean <- dat %>%
  dplyr::select(-s1_start_centre) %>%
  mutate(pp_prior, pp_after) %>%
  filter(pp_prior > 5)

saveRDS(dat_clean, "dat/plots_pp.rds")

ggplot() + 
  geom_point(data = dat_clean, aes(x = pp_prior, y = s1_start))

ggplot() + 
  geom_point(data = dat_clean, aes(x = pp_after, y = s1_end))
