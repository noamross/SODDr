#' Given a lattice and a matched set of plots with data, populate the rest of
#' lattice 
#' @import geoRglm
PopulateLattice <- function(lattice, plotdata=NULL, method="poisson") {
  
  if(method=="poisson") {

  # Calculate the means of each class
    
  # Assign the matched lattice points values from the data set
    
  # Calculate 
    
  }
  
  if(method="krige") {
    
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