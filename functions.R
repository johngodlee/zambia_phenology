# Define variable name lookup
pred_lookup <- c(
  eff_rich = "Diversity",
  evenness = "Evenness",
  n_stems_ge10_ha = "Stem density",
  cum_precip_seas = "Wet season precip.",
  cum_precip_pre = "Pre-green-up precip.",
  cum_precip_end = "Pre-senescence precip.",
  cluster = "Vegetation type",
  diurnal_temp_range = "Diurnal dT",
  diam_quad_mean = "Stem diameter",
  diam_cov = "Stem diameter variance",
  Detarioideae = "Detarioid BA")

resp_lookup <- c(
  cum_vi = "Cumulative EVI",
  s1_length = "Season length",
  s1_green_rate = "Green-up rate",
  s1_senes_rate = "Senescence rate",
  start_lag = "Green-up lag",
  end_lag = "Senescence lag")

# Define cluster name lookup
clust_lookup <- c(
  `1` = "1",
  `2` = "2",
  `3` = "3",
  `4` = "4",
  `5` = "5",
  `6` = "6")

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

#' Generate a valid UTM WGS84 proj4string given a UTM zone
#'
#' @param x character vector defining UTM zones
#'
#' @return UTM proj4string character vector
#' 
UTMProj4 <- function(x){
  unlist(lapply(1:length(x), function(y) {
    paste0(
      "+proj=utm +zone=",
      gsub("[A-z]", "", as.character(x[y])),
      ifelse(gsub("[0-9]", "", as.character(x[y])) == "S", " +south", ""),
      " +ellps=WGS84")
  }))
}

#' Format numbers for LaTeX
#'
#' @param x numeric atomic
#' @param digits integer number of digits to round to
#' @param method method of rounding either "round" or "signif"
#'
#' @return formatted character strign
#' 
numFormat <- function(x, digits = 2, method = "round"){
  sprintf(paste0("%.",digits,"f"),
    if (method == "round") {
      round(x, digits = digits)
    } 
    else if (method == "signif") {
      signif(x, digits = digits)
    })
}

#' Format a p-value for use in LaTeX
#'
#' @param p atomic vector p-value
#' @param digits number of decimal places
#'
#' @return character string
#' 
#' @examples
#' p <- 0.04592
#' pFormat(p)
#' 
pFormat <- function(p, print_p = TRUE, digits = 2){
  unlist(lapply(p, function(x) {
    if (x < 0.01) {
      if (print_p) {
        return("p<0.01")
      } else {
        return("<0.01")
      }
    } 
    else if (x < 0.05) {
      if (print_p) {
        return("p<0.05")
      } else {
        return("<0.05")
      }
    }
    else {
      out <- as.character(numFormat(x, digits = digits))
      if (print_p) {
        return(paste0("p = ", out))
      } else {
        return(out)
      }
    }
  }))
}

#' Create LaTeX newcommand output
#'
#' @param x atomic vector to export
#' @param name LaTeX variable name 
#'
#' @return string
#' 
commandOutput <- function(x, name){ 
  paste0("\\newcommand{\\",
    ifelse(missing(name), deparse(substitute(x)), name), 
    "}{", 
    x, 
    "}"
  )
}

corrPlot <- function(x, col = c("blue", "white", "red"), ...) {
  corr <- psych::corr.test(x, ...) 
  corr_ci <- data.frame(raw.lower = corr$ci$lower, raw.r = corr$ci$r, 
    raw.upper = corr$ci$upper, raw.p = corr$ci$p, 
    adj.lower = corr$ci.adj$lower.adj, adj.upper = corr$ci.adj$upper.adj)
  corr_ci$vars <- row.names(corr_ci)
  corr_ci$conf_x <- unlist(sapply(1:(length(x)-1), function(i){
      c(1:(length(x)-1))[i:(length(x)-1)]
    })) + 1
  rev_mat <- (length(x)-1):1
  corr_ci$conf_y <- unlist(sapply(1:(length(x)-1), function(i){
      rep(i, times = rev_mat[i])
    }))
  n_seq <- 2:length(x)
  corr_ci$y_var <- unlist(sapply(1:(length(x)-1), function(i){
      rep(row.names(corr[[1]])[i], rev_mat[i])
    }))
  corr_ci$x_var <- unlist(sapply(1:(length(x)-1), function(i){
      row.names(corr[[1]])[n_seq[i]:length(x)]
    }))
  corr_ci$x_var <- factor(corr_ci$x_var, levels = unique(corr_ci$x_var))
  corr_ci$y_var <- factor(corr_ci$y_var, levels = unique(corr_ci$y_var))
  corr_ci$conf <- (corr_ci$raw.lower > 0) == (corr_ci$raw.upper > 0)
  corr_ci$raw.r <- round(corr_ci$raw.r, 2)

  ggplot2::ggplot() + 
    ggplot2::geom_tile(data = corr_ci, 
      ggplot2::aes(x = x_var, y = y_var, 
        fill = raw.r), colour = "black") + 
    ggplot2::geom_text(data = corr_ci, 
      ggplot2::aes(x = x_var, y = y_var, label = raw.r),
      size = 3) + 
    ggplot2::geom_point(data = corr_ci[corr_ci$conf == FALSE,], 
      ggplot2::aes(x = x_var, y = y_var), 
      fill = NA, colour = "black", shape = 21, size = 11) + 
    ggplot2::scale_fill_gradient2(name = "r", 
      low = col[1], mid = col[2], high = col[3]) + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x = "", y = "") + 
    ggplot2::coord_equal() +
    ggplot2::theme(legend.position = "none")
}

# Get convex hull of 2D points
findHull <- function(dat, x, y, group = NULL) {
  if (is.null(group)) {
    out <- dat[chull(dat[[x]], dat[[y]]), c(x, y)]

  } else {
    hulls <- by(dat, dat[[group]], function(i) {
      i[chull(i[[x]], i[[y]]), c(x, y, group)]
    })

    out <- do.call(rbind, hulls)
  }

  return(out)
}


# Decompose time series with overlap
seasonGet <- function(x, min_date, max_date, date = "date") {
  # Format date as date
  x[[date]] <- as.Date(x[[date]])

  # Subset data between max and min dates
  out <- x[x[[date]] > as.Date(min_date) & x[[date]] < as.Date(max_date),]

  if (nrow(out) > 0) {
    # Get the "season" of the data, i.e. the year the measurements start
    out$season <- min(format(out[[date]], "%Y"))

    # Get the days from the min date
    out$days <- as.numeric(out[[date]] - as.Date(min_date))

    # Get year
    out$year <- format(out[[date]], "%Y")

    # Recenter days on start of year in middle of growing season
    out$doy <- as.numeric(out$date - as.Date(paste0(format(as.Date(max_date), "%Y"), "-01-01")))
    
    return(out)
  }
}

#' Generate a species by site abundance matrix
#'
#' @param x dataframe of individual records
#' @param site_id column name string of site IDs
#' @param species_id column name string of species names
#' @param fpc optional column name string of sampling weights of each record, 
#'     between 0 and 1 
#' @param abundance optional column name string with an alternative abundance 
#'     measure such as biomass, canopy cover, body length
#'
#' @return dataframe of species abundances (columns) per site (rows)
#' 
#' @examples
#' x <- data.frame(site_id = rep(c("A", "B", "C"), each = 3), 
#'   species_id = sample(c("a", "b", "c", "d"), 9, replace = TRUE), 
#'   fpc = rep(c(0.5, 0.6, 1), each = 3), 
#'   abundance = seq(1:9))
#' abMat(x, "site_id", "species_id")
#' abMat(x, "site_id", "species_id", "fpc")
#' abMat(x, "site_id", "species_id", "fpc", "abundance")
#' 
#' @export
#' 
abMat <- function(x, site_id, species_id, fpc = NULL, abundance = NULL) {
  # If no fpc or abundance, make 1
  if (is.null(fpc)) {
    x$fpc <- 1
  } else {
  	x$fpc <- x[[fpc]]
  }
  if (is.null(abundance)) {
    x$abundance <- 1 
  } else {
  	x$abundance <- x[[abundance]]
  }

  x$abundance <- x$abundance / x$fpc

  # Count number of species and sites
  comm_df <- aggregate(x$abundance, by = list(x[[site_id]], x[[species_id]]), 
    simplify = FALSE, drop = FALSE, FUN = sum)

  # Replace NULL with zero
  comm_df$x <- unlist(lapply(comm_df$x, function(y) {
      if(is.null(y)) {
        0
      } else {
        y
      }
    }))
  
  # Make names tidy
  names(comm_df) <- c("x","y","z")
  comm_df$x <- factor(comm_df$x)
  comm_df$y <- factor(comm_df$y)

  # Spread to matrix
  comm <- with(comm_df, {
    out <- matrix(nrow = nlevels(x), ncol = nlevels(y),
      dimnames = list(levels(x), levels(y)))
    out[cbind(x, y)] <- z
    out
  })

  comm <- as.data.frame(comm)

  return(comm)
}


#' Calculate the quadratic mean
#'
#' @param x vector 
#'
#' @return quadratic mean of vector
#' 
quadMean <- function(x, ...) {
  sqrt(mean(x^2, ...))
}
