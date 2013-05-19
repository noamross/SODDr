#' Create a hexagonal grid
#' 
#' Given a data frame or a SpatialPointsDataFrame, return a hexagonal grid
#' of points covering its extent
#' 
#' @param data A data frame or SpatialPointsDataFrame containing coordinates
#' data
#' @param hex.area the area of each hexagonal tile
#' @param xy A two-element vector with the names of the columns in the data
#' frame indicating coordinates.  Not needed for a SpatialPointsDataFrame
#' @import sp
#' @export
hex.grid <- function(data, hex.area, xy=c("Easting", "Northing")) {
    if (class(data) != "SpatialPointsDataFrame") {
      coordinates(data) <- xy
    }
    hex.pts <- spsample(data, type = "hexagonal", 
                     cellsize = sqrt(hex.area*2/(sqrt(3))), 
                     bb=bbox(data))
    locations <- coordinates(hex.pts)
    rownames(locations) <- 1:nrow(locations)
    colnames(locations) <- colnames(coordinates(data))
    return(locations)
}