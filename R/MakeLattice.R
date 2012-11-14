#'Generate a lattice of equally spaced locations and coordinates
#'@export
MakeLattice <- function(nx, ny, dist=1) {
  locations <- cbind(location = 1:(nx*ny), 
                     x = rep(seq(from=0, by=dist, length.out=nx), each=ny),
                     y = rep(seq(from=0, by=dist, length.out=ny), times=nx)
  )
  return(locations)
}