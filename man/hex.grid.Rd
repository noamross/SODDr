\name{hex.grid}
\alias{hex.grid}
\title{Create a hexagonal grid}
\usage{
  hex.grid(data, hex.area, xy = c("Easting", "Northing"))
}
\arguments{
  \item{data}{A data frame or SpatialPointsDataFrame
  containing coordinates data}

  \item{hex.area}{the area of each hexagonal tile}

  \item{xy}{A two-element vector with the names of the
  columns in the data frame indicating coordinates.  Not
  needed for a SpatialPointsDataFrame}
}
\description{
  Given a data frame or a SpatialPointsDataFrame, return a
  hexagonal grid of points covering its extent
}

