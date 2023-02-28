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


