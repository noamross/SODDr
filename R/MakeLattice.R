#'Generate a lattice of equally spaced locations and coordinates
#'@export
MakeLattice <- function(nx, ny, dist=1) {
  locations <- cbind(location = 1:(nx*ny), 
                     x = rep(seq(from=0, by=dist, length.out=nx), each=ny),
                     y = rep(seq(from=0, by=dist, length.out=ny), times=nx)
  )
  class(locations) <- c(class(locations), "lattice")
  attr(locations, "gridsize") <- dist
  return(locations)
}

#'Generate a lattice of points equally spaced in the centers of a hexagonal
#'lattice
#'@export
MakeHexLattice <- function(nx,ny,dist=1) {
  locations <- cbind(location = 1:(nx*ny),
                     x = sort(c(rep(seq(from=0, by=dist, length.out=nx),
                                    each=ceiling(ny/2)),
                                rep(seq(from=dist/2, by=dist, length.out=nx),
                                    each=floor(ny/2)))),
                     y = rep(c(seq(from=0, by = dist*sqrt(3), 
                                   length.out=ceiling(ny/2)),
                               seq(from=dist*sqrt(3)/2, by=dist*sqrt(3),
                                   length.out=floor(ny/2))),
                             times=nx))
  class(locations) <- c(class(locations), "lattice")
  attr(locations, "gridsize") <- dist
  return(locations)
}

