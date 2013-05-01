#' Given a lattice and a matched set of plots with data, populate the rest of
#' lattice
#' @param lattice
#' @param plotdata A data frame listing plot ids and counts of each class 
#' @import geoRglm
PopulateLattice <- function(lattice, plotdata=NULL, method="poisson") {
  
  if(method == "poisson") {

  # Calculate the means of each class
    class.means <- colMeans(plotdata[2:ncols(plotdata)])
  # Assign the matched lattice points values from the data set
    
  # Calculate 
    
  }
  
  if(method=="krige") {
    
  # Fit a krige spatial model to the species distributions
  krige.model <- pois.krige(coords=lattice[locations$lattice.nn,2:3],
                            data=locations[TODO],
                            locations=lattice[-locations$lattice.nn,2:3],
                            mcmc.input=,
                            krige=list(type.krige="ok",
                            output=list(sim.posterior=FALSE, sim.predict=TRUE
                                        messages=FALSE))
  # Simulate from the model for each point on the lattice
    
  }
}