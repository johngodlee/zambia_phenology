# Generate abundance matrices for sites
# John L. Godlee (johngodlee@gmail.com)
# Last updated: 2023-02-28

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

# Import data
trees <- readRDS("./dat/trees.rds")

# Create tree species abundance matrices by plot cluster
ab_plot_mat <- abMat(trees, 
  site_id = "plot_id", 
  species_id = "species_name_clean") 
ab_plot_mat <- ab_plot_mat[,names(ab_plot_mat) != "Indet indet"]

ba_clust_mat <- abMat(trees, 
  site_id = "plot_cluster", 
  species_id = "species_name_clean", 
  abundance = "ba")
ba_clust_mat <- ba_clust_mat[,names(ba_clust_mat) != "Indet indet"]

# Write files
saveRDS(ba_clust_mat, "dat/ba_clust_mat.rds")
saveRDS(ab_plot_mat, "dat/ab_plot_mat.rds")

