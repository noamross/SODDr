\name{InitModel}
\alias{InitModel}
\title{Generate a Model from which to draw a random forest based on plot data}
\usage{
  InitModel(data, plot.area = 500,
    model.type = c("glm", "gam"), Count = "Count",
    form = NULL, xy = c("Easting", "Northing"),
    p4string = "+proj=utm +zone=10 +ellps=intl +units=m")
}
\description{
  Generate a Model from which to draw a random forest based
  on plot data
}

