## Performance flags

parallel.flag <- FALSE
progress.flag <- !parallel.flag && interactive()

# 1.  Set up species/classes
#' Import a CSV of parameters, process, and return objects to parent environment fadskjhf asjdhgf adgjhadsk
#' @import gdata
  treeparms.df <- read.csv("test_tree_parms.csv", stringsAsFactors=FALSE)
  classes <- 1:nrow(treeparms.df)
    names(classes) <- paste(treeparms.df$species, ",", treeparms.df$ageclass,sep="")
  treeparms.df <- cbind(class=classes, treeparms.df)
  n.species <- length(unique(treeparms.df$species))
  n.classes <- nrow(treeparms.df)
  classes <- 1:n.classes
  ageclasses <-as.vector(table(treeparms.df$species))
    names(ageclasses) <- 1:n.species
  waifw <- as.matrix(treeparms.df[,matchcols(treeparms.df, with="waifw[0-9]+"),])
    dimnames(waifw) <- list(names(classes),names(classes))
  recruit.vec <- as.vector(rbind(treeparms.df$S.recruit, treeparms.df$I.recruit))
  mort.vec <- as.vector(rbind(treeparms.df$S.mortality, treeparms.df$I.mortality))
  trans.vec <- as.vector(rbind(treeparms.df$S.transition, treeparms.df$I.transition))
  resprout.vec <- as.vector(rbind(treeparms.df$S.resprout, treeparms.df$I.resprout))
  recover <- treeparms.df$recover
# 2. Set up locations.  In this case, a lattice

#'Generate a lattice of equally spaced locations and coordinates
MakeLattice <- function(nx, ny, dist=1) {
  locations <- cbind(location = 1:(nx*ny), 
                     x = rep(seq(from=0, by=dist, length.out=nx), each=ny),
                     y = rep(seq(from=0, by=dist, length.out=ny), times=nx)
                )
  return(locations)
}

locations <- MakeLattice(5,5,1)
n.locations <- nrow(locations)
  
#'Generate a distance matrix from locations.

distance.matrix <- as.matrix(dist(locations[,c("x","y")]))
require(plyr)
  
# Generate dispersal matrices for each class
spread.matrices <- daply(treeparms.df,"class", function(x) {
  do.call(x$kernel.fn, 
          unname(c(list(distance.matrix), 
                   as.list(na.omit(x[matchcols(x, "kernel.par[0-9]+")]))
          ))
  )
})
  

# Set simulation parameters

time.steps <- 1:50

# Create transition matrix
# TODO: Make this into a function
  
tran.mat <- matrix(0, nrow=n.classes*2, ncol=n.classes*2)
diags <- row(tran.mat) - col(tran.mat)
fec.mat <- tran.mat
for(i in 1:n.species) {
  classindex <- 1:(2*ageclasses[i]) + (sum(ageclasses[0:(i-1)])*2)
  fec.mat[classindex[1],classindex] <- (recruit.vec + resprout.vec*mort.vec)[classindex]  #Fecundities + Death X Resprout probabilties
}
diag(tran.mat) <- 1 - trans.vec - mort.vec
  
tran.mat[diags==1] <- diag(tran.mat)[-(nrow(tran.mat))] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
tran.mat[diags==-1] <- diag(tran.mat)[-1] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]

tran.mat[diags==2] <- trans.vec[1:(n.classes*2 - 2)]
tran.mat[diags==3] <- trans.vec[1:(n.classes*2 -3)] * rep(c(1,0),n.classes)[1:(n.classes*2 - 3)]
tran.mat[diags==1] <- tran.mat[diags==1] + trans.vec[1:(n.classes*2-1)] * rep(c(0,1), n.classes)[1:(n.classes*2 - 1)]
  
# Step update is (trans.mat * force.mat + fec.mat) * population
  
# Initiatte
  
pop <- array(NA, dim=c(n.locations, 2*n.classes, length(time.steps)))   #TODO: Standardize use of n.classes
pop.init <- matrix(data=c(0.125,0), nrow=n.locations, ncol=(2*n.classes), byrow=TRUE)
pop[,,1] <- pop.init

spore.burden <- matrix(NA, nrow=n.classes, ncol=n.locations)
for(time in time.steps[-(length(time.steps))]) {
  
# First act in simulation step.  Given population at each location, calculate spore burden at each location
  
  for(class in classes) {
      spore.burden[class,] <- pop[,class*2,time] %*% spread.matrices[class,,]
  }
  
  force.infection <- waifw %*% spore.burden
  real.recovery <- (1-force.infection) * recover

  for(location in 1:n.locations) {
    force <- force.infection[,location]
    real.rec <- real.recovery[,location]
    force.matrix <- kronecker(rep(1,n.classes),
                              rbind(as.vector(rbind(1-force,real.rec)),
                              as.vector(rbind(force, 1-real.rec))))
    trans.mat <- tran.mat*force.matrix + fec.mat
    pop[location,,time + 1] <- trans.mat %*% pop[location,,time]
  }
}
  
  # Calculate spore burden matrix
  

  


# spore.burden is a matrix at time t of the spores in each 
# cell coming FROM trees of each class in all other cells

spore.burden <- matrix(NA, prod(lattice.dims), sum(size.classes))

# To fill, multiply infected population of each class by weighted.matrix
infected.pop <- pop.array[,seq(2,classes*2,2),t]
spore.burden <- infected.pop %*% weighted.matrix
# This can be multiplied by waifm to get the force of infection on each species at each site

lambda <- spore.burden * waifw


## Use popbio::multiresultm for demographic stochasticity

rowSums
