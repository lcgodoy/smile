#' Find phi parameter for the Exponential spatial auto-correlation function
#'
#' Function designed to find the phi paramter such that the correlation
#' between points wihtin a given distance \code{d} is at most a given
#' value.
#'
#' @param d maximun distance for spatial dependence
#' @param kappa parameter
#' @param range ? Marcos help needed here
#' @param cut ? Marcos help needed here
#'
#' @return \code{real number}
#' @export
#'
find_phi <- function(d, kappa = 0.5, 
                     range = c(1e-04, 1), 
                     cut = 0.05) {    
  out <- stats::uniroot(f = function(x, d, kappa, cut) {      
    out <- 1/(2^(kappa - 1)*gamma(kappa)) * (d/x)^kappa * besselK(d/x,kappa)
    out <- out - cut    
  },
  interval = range, d = d, kappa = kappa, cut = cut)    
  out$root 
}
