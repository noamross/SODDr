
DensityDependence.Rfun <- function(pop, compete, K) { 
  (K - sum(pop*compete))/K
}


require(compiler)
#' Produced multiplier for new recruits given population and competition vectors
#'@import compiler
#'@export
DensityDependence <- cmpfun(DensityDependence.Rfun)
