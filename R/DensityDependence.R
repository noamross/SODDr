
DensityDependence.Rfun <- function(pop, compete) { 
  1 - sum(pop*compete)
}


require(compiler)
#' Produced multiplier for new recruits given population and competition vectors
#'@import compiler
#'@export
DensityDependence <- cmpfun(DensityDependence.Rfun)
