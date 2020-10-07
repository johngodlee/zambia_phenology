# Theme colours
pal <- c("lightseagreen", "#DE6400", "dodgerblue", "tomato", "grey", "#E0E0E0")

# Cluster colours
clust_pal <- c("#5cb930", "#6e49cc", "#49b9da", "#da2c57")

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

#' Get valid UTM zone from latitude and longitude in WGS84 decimal degrees
#'
#' @param x vector of longitude coordinates in decimal degrees
#' @param y vector of latitude coordinate in decimal degrees
#'
#' @return Vector of UTM zones for each latitude-longitude pair
#' 
#' @export
#' 
latLong2UTM <- function(x, y) {
  unlist(lapply(1:length(x), function(z) {
    paste((floor((x[z] + 180) / 6) %% 60) + 1,
      ifelse(y[z] < 0, "S", "N"),
      sep = "")
  }))
}

#' Generate a valid UTM WGS84 proj4string given a UTM zone
#'
#' @param x character vector defining UTM zones
#'
#' @return UTM proj4string character vector
#' 
#' @export
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

#' 
#'
#' @param sort logical, if true rows are sorted according to vector length on first two axes
#'
#' @return 
#' 
#' @examples
#' 
#' 
#' @export
#' 
pcoaArrows <- function(given_pcoa, trait_df, sort = FALSE) {
    n <- nrow(trait_df)
    points.stand <- scale(given_pcoa$vectors)
    
    # Compute covariance of variables with all axes
    S <- cov(trait_df, points.stand)
    
    # Select only positive eigenvalues
    pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
    
    # Standardize value of covariance (see Legendre & Legendre 1998)
    U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
    colnames(U) <- colnames(given_pcoa$vectors)
    
    # Add values of covariances inside object
    given_pcoa$U <- U

    if (sort) {
      # Compute arrow vector lengths
      dists <- list()
      for (i in 1:nrow(given_pcoa$U)) {
        dists[i] <- as.vector(dist(t(matrix(c(given_pcoa$U[i,1:2], 0,0), 2))))
      }

      # Order arrows dataframe by vector length 
      given_pcoa$U <- as.data.frame(given_pcoa$U[order(unlist(dists), 
          decreasing = TRUE),])
    }
    
    return(given_pcoa)
}

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


commandOutput <- function(x, name){ 
  paste0("\\newcommand{\\",
    ifelse(missing(name), deparse(substitute(x)), name), 
    "}{", 
    x, 
    "}"
  )
}

#' spaMM model effect sizes 
#'
#' @param mod spaMM model object
#' @param intercept logical, should the intercept term be included?
#'
#' @return ggplot2 object
#' 
spammEff <- function(mod, intercept = FALSE) {
  capture.output(mod_summ <- summary(mod))
  slopes <- data.frame(var = row.names(mod_summ$beta_table), 
    est = mod_summ$beta_table[,1], 
    se = mod_summ$beta_table[,2])
  row.names(slopes) <- NULL
  if (!intercept) {
    slopes <- slopes[slopes$var != "(Intercept)",]
  }
  slopes$var <- factor(slopes$var, levels = slopes$var[order(nrow(slopes):1)])
  return(slopes)
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

#' Extract predicted grid from spaMM model over range of observations
#'
#' @param fitobject spaMM model object
#' @param gridSteps number of levels of xy grid
#'
#' @return dataframe of xy grid with predicted values
#' 
spammMapExtract <- function (fitobject, gridSteps = 200) {
  # Obtain coordinates
  info_olduniqueGeo <- attr(fitobject, "info.uniqueGeo")
  coordinates <- unique(unlist(lapply(info_olduniqueGeo, colnames)))

  pred <- predict(fitobject, binding = "fitted")
  form <- spaMM::formula.HLfit(fitobject, which = "hyper")
  map.formula <- as.formula(paste(attr(pred, "fittedName"),
    " ~ 1 + ", paste(spaMM:::.findSpatial(form), collapse = " + ")))
  smoothObject <- fitme(map.formula, data = pred, 
    fixed = list(phi = 1e-05), method = "REML")
  smoo <- predict(smoothObject, binding = "dummy")
  x <- smoo[, coordinates[1]]
  y <- smoo[, coordinates[2]]

  marg <- 1/20
  xrange <- range(x)
  margex <- (xrange[2] - xrange[1]) * marg
  xrange <- xrange + margex * c(-1, 1)

  yrange <- range(y)
  margey <- (yrange[2] - yrange[1]) * marg
  yrange <- yrange + margey * c(-1, 1)

  xGrid <- seq(xrange[1], xrange[2], length.out = gridSteps)
  yGrid <- seq(yrange[1], yrange[2], length.out = gridSteps)
  newdata <- expand.grid(xGrid, yGrid)
  colnames(newdata) <- coordinates

  newdata$val <- predict(smoothObject, newdata = newdata, 
    variances = list(), control = list(fix_predVar = FALSE))
  return(newdata)
}


#' Calculate Area Under Curve (AUC) using trapesium approximation
#'
#' @param x vector of x values
#' @param y vector of y values
#'
calcAUC <- function(x, y) {
  sum( diff(x) * (head(y,-1) + tail(y,-1)) ) / 2
}

# Run many iterations of NMDS to check optimal number of dimensions
NMDS.scree <- function(x, dims = 10, ...) {
  # Create dimensions vector
  if (length(dims) == 1) {
    dim_vec <- seq(dims) 
  } else {
    dim_vec <- dims
  }

  # Create list of metaMDS objects
  meta_list <- lapply(dim_vec, function(y) {
    metaMDS(x, k = y, ...)
  })

  # Extract stress values
  stress_vec <- unname(unlist(lapply(meta_list, `[`, "stress")))

  # Create plot
  plot(dim_vec, stress_vec,
    xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
}

# Build ellipses from vegan::ordiellipse
covEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- seq(0, npoints) * 2 * pi / npoints
  circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(circle %*% chol(cov)))
}

# Here's the stat_
StatBag <- ggproto("Statbag", Stat,
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
GeomBag <- ggproto("GeomBag", Geom,
                   draw_group = function(data, panel_scales, coord) {
                     n <- nrow(data)
                     if (n == 1) return(zeroGrob())

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

                   default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1,
                                     alpha = NA, prop = 0.5),

                   handle_na = function(data, params) {
                     data
                   },

                   required_aes = c("x", "y"),

                   draw_key = draw_key_polygon
)
