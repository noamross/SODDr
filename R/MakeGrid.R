#'Generate a lattice of equally spaced locations and coordinates
#'@export
MakeGrid <- function(nx=NULL, ny=NULL, dist=1, type="square", area=NULL,
                     origin=c(0,0), around=NULL, xy=c("Easting", "Northing")) {
                       
                       
  if(!is.null(around)) {
    origin <- c(min(around[xy[1]]), min(around[xy[2]]))
    maxdim <- c(max(around[xy[1]]), max(around[xy[2]]))
  }
  
  if(!is.null(area) & type=="square") 
  if(!is.null(area) & type=="hex") dist <- sqrt(area*2/sqrt(3))
  
  if(type=="square") {
    if(!is.null(area)) {
      dist <- sqrt(area)
      nx <- floor((maxdim-origin)[1]/dist)
      ny <- floor((maxdim-origin)[2]/dist)
    }
    x <- rep(seq(from=0, by=dist, length.out=nx), each=ny) + origin[1]
    y <- rep(seq(from=0, by=dist, length.out=ny), times=nx) + origin[2]
  } else if(type=="hex") {
    if(!is.null(area)) {
      dist <- sqrt(area*2/sqrt(3))
      nx <- floor((maxdim-origin)[1]/dist)
      ny <- floor((maxdim-origin)[2]/(dist*sqrt(3)/2))
    }
    x <- sort(c(rep(seq(from=0, by=dist, length.out=nx), each=ceiling(ny/2)),
               rep(seq(from=dist/2, by=dist, length.out=nx), each=floor(ny/2))
               )) + origin[1]
    y <- rep(c(seq(from=0, by = dist*sqrt(3), length.out=ceiling(ny/2)),
              seq(from=dist*sqrt(3)/2, by=dist*sqrt(3), length.out=floor(ny/2))
              ) + origin[2], times=nx)
  }
    
  locations <- cbind(x, y)
  dimnames(locations) <- list(Location=1:nrow(locations),
                              Coordinates=xy)
  class(locations) <- c(class(locations), paste0(type,"_grid"))
  attr(locations, "gridsize") <- dist
  if(!is.null(area)) attr(locations, "cell_area") <- area
  return(locations)
}