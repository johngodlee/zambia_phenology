# Define variable name lookup
pred_lookup <- c("Richness", "Evenness", "Stem density",  
  "MAP", "Wet season precip.", "Pre-green-up precip.", "Pre-senescence precip.",
    "Vegetation type", "Diurnal dT")
names(pred_lookup) <- c("eff_rich", "evenness", "n_stems_gt10_ha",
  "map", "cum_precip_seas", "cum_precip_pre", "cum_precip_end",
  "cluster", "diurnal_temp_range")

resp_lookup <- c("Cumulative EVI", "Season length", 
  "Green-up rate", "Senescence rate",
  "Green-up lag", "Senescence lag")
names(resp_lookup) <- c("cum_vi", "s1_length", 
  "s1_green_rate", "s1_senes_rate", 
  "start_lag", "end_lag")

# Define cluster name lookup
clust_lookup <- c("1", "2", "3", "4", "5", "6")
names(clust_lookup) <- c("1", "2", "3", "4", "5", "6")

# Theme colours
pal <- c("lightseagreen", "#DE6400", "dodgerblue", "tomato", "grey", "#E0E0E0")

# Cluster colours
clust_pal <- c("lightseagreen", "#1b9e77","#d95f02","#7570b3","#e7298a",
  "#66a61e","#49b9da","#a6761d", "tomato", "grey", "navy", "forestgreen", 
  "darkgoldenrod", "black")

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

# ggplot2 convex hulls
# Here's the stat_
StatBag <- ggplot2::ggproto("Statbag", ggplot2::Stat,
                   compute_group = function(data, scales, prop = 0.5) {

                     #################################
                     #################################
                     # originally from aplpack package, plotting functions removed
                     plothulls_ <- function(x, y, fraction, n.hull = 1,
                                            col.hull, lty.hull, lwd.hull, density=0, ...){
                       # function for data peeling:
                       # x,y : data
                       # fraction.in.inner.hull : max percentage of points within the hull to be drawn
                       # n.hull : number of hulls to be plotted (if there is no fractiion argument)
                       # col.hull, lty.hull, lwd.hull : style of hull line
                       # plotting bits have been removed, BM 160321
                       # pw 130524
                       if(ncol(x) == 2){ y <- x[,2]; x <- x[,1] }
                       n <- length(x)
                       if(!missing(fraction)) { # find special hull
                         n.hull <- 1
                         if(missing(col.hull)) col.hull <- 1
                         if(missing(lty.hull)) lty.hull <- 1
                         if(missing(lwd.hull)) lwd.hull <- 1
                         x.old <- x; y.old <- y
                         idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
                         for( i in 1:(length(x)/3)){
                           x <- x[-idx]; y <- y[-idx]
                           if( (length(x)/n) < fraction ){
                             return(cbind(x.hull,y.hull))
                           }
                           idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx];
                         }
                       }
                       if(missing(col.hull)) col.hull <- 1:n.hull
                       if(length(col.hull)) col.hull <- rep(col.hull,n.hull)
                       if(missing(lty.hull)) lty.hull <- 1:n.hull
                       if(length(lty.hull)) lty.hull <- rep(lty.hull,n.hull)
                       if(missing(lwd.hull)) lwd.hull <- 1
                       if(length(lwd.hull)) lwd.hull <- rep(lwd.hull,n.hull)
                       result <- NULL
                       for( i in 1:n.hull){
                         idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
                         result <- c(result, list( cbind(x.hull,y.hull) ))
                         x <- x[-idx]; y <- y[-idx]
                         if(0 == length(x)) return(result)
                       }
                       result
                     } # end of definition of plothulls
                     #################################


                     # prepare data to go into function below
                     the_matrix <- matrix(data = c(data$x, data$y), ncol = 2)

                     # get data out of function as df with names
                     setNames(data.frame(plothulls_(the_matrix, fraction = prop)), nm = c("x", "y"))
                     # how can we get the hull and loop vertices passed on also?
                   },

                   required_aes = c("x", "y")
)

# Here's the stat_ function
#' @inheritParams ggplot2::stat_identity
#' @param prop Proportion of all the points to be included in the bag (default is 0.5)
stat_bag <- function(mapping = NULL, data = NULL, geom = "polygon",
                     position = "identity", na.rm = FALSE, show.legend = NA, 
                     inherit.aes = TRUE, prop = 0.5, alpha = 0.3, ...) {
  layer(
    stat = StatBag, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...)
  )
}

# here's the geom_
geom_bag <- function(mapping = NULL, data = NULL,
  stat = "identity", position = "identity",
  prop = 0.5, 
  alpha = 0.3,
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = StatBag,
      geom = GeomBag,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm,
        alpha = alpha,
        prop = prop,
        ...
      )
    )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomBag <- ggplot2::ggproto("GeomBag", ggplot2::Geom,
draw_group = function(data, panel_scales, coord) {
  n <- nrow(data)
  if (n == 1){ return(zeroGrob()) }

  munched <- coord_munch(coord, data, panel_scales)

  # Sort by group to make sure that colors, fill, etc. come in same order
  munched <- munched[order(munched$group), ]

  # For gpar(), there is one entry per polygon (not one entry per point).
  # We'll pull the first value from each group, and assume all these values
  # are the same within each group.
  first_idx <- !duplicated(munched$group)
  first_rows <- munched[first_idx, ]

  ggplot2:::ggname("geom_bag",
    grid:::polygonGrob(munched$x, munched$y, default.units = "native",
      id = munched$group,
      gp = grid::gpar(
        col = first_rows$colour,
        fill = alpha(first_rows$fill, first_rows$alpha),
        lwd = first_rows$size * .pt,
        lty = first_rows$linetype
      )
    )
  )
  },
  default_aes = ggplot2::aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1,
    alpha = NA, prop = 0.5),

  handle_na = function(data, params) {
    data
  },

  required_aes = c("x", "y"),

  draw_key = ggplot2::draw_key_polygon
)

# Decompose time series with overlap
seasonGet <- function(x, min_date, max_date, date = "date") {
  # Format date as date
  x[[date]] <- as.Date(x[[date]])

  # Subset data between max and min dates
  out <- x[x[[date]] > as.Date(min_date) & x[[date]] < as.Date(max_date),]

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
