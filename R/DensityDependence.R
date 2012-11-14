
DensityDependence.Rfun <- function(pop, space) { 
  1 - sum(pop*space)
}


require(compiler)
#' Produced multiplier for new recruits given population and space vectors
#'@import compiler
#'@export
DensityDependence <- cmpfun(DensityDependence.Rfun)
