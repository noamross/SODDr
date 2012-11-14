#'@import compiler
DensityDependence.Rfun <- function(pop, space) { 
  1 - sum(pop*space)
}

#' Produced multiplier for new recruits given population and space vectors
DensityDependence <- cmpfun(DensityDependence.Rfun)
