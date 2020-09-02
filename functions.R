
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

commandOutput <- function(x, name){ 
  paste0("\\newcommand{\\",
    ifelse(missing(name), deparse(substitute(x)), name), 
    "}{", 
    x, 
    "}"
  )
}

#' Interval plot of spaMM model effect sizes 
#'
#' @param mod spaMM model object
#' @param intercept logical, should the intercept term be included?
#'
#' @return ggplot2 object
#' 
spammEffPlot <- function(mod, intercept = FALSE) {
  capture.output(mod_summ <- summary(mod))
  slopes <- data.frame(var = row.names(mod_summ$beta_table), 
    est = mod_summ$beta_table[,1], 
    se = mod_summ$beta_table[,2])
  if (!intercept) {
    slopes <- slopes[slopes$var != "(Intercept)",]
  }
  out <- ggplot() + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point(data = slopes, aes(x = est, y = var)) + 
    geom_errorbarh(data = slopes, 
      aes(xmin = est - se, xmax = est + se, y = var), height = 0.1) + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line(colour = "#E0E0E0"),
      axis.text = element_text(size = 12)) + 
    labs(x = "Effect slope", y = "Factor")

  print(out)
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
