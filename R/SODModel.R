# TODO: Check that input dimensions match up
#  - Locations and classes against initial conditions
#  - Create separate function error_check
  
# Make robust to single location, single class, single replicate
# First check if dimnames are right, and use aperm
# Otherwise use order of array
# Make sure dimnames flow through model, 
# Make output higher-dimensional array?

# Make this interruptible? How, saving results to disk each time?
  
# Make base function take vector of parameters, with matrix for replicate runs
# Check dimensions of parameters and initial conditions, if same do that number of runs
# If different then multiply.

# Add run time attribute

# Use this for progress bar: http://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply


#TODO: Pre-calculate dispersal matrix for all runs (and other overhead, e.g.
#turning parameters into a vector)

#'Run the SOD model
#'@param parallel Run the simulation in parallel or not?
#'@param ... arguments to be passed to SODModel
#'@import foreach doRNG iterators abind plyr
#'@export
SODModel <- function(parms, times, locations, init, reps=1, lambda.ex = 0, K=50,
                     stochastic.e=NULL, stochastic.d=FALSE, 
                     verbose=interactive(), parallel=FALSE) {
  
  init <- abind(rep(list(init),each=reps), along=3)
  n.runs <- dim(init)[3]

  if(verbose) {
    message(n.runs, " Simulations of ",
            tail(times, 1) - times[1] + 1, " time steps each")
    if(!parallel) {
        z <- plyr:::txtTimerBar((tail(times, 1) - times[1] + 1)*n.runs)
    }
  }
  
  `%op%` <- if (parallel) `%dopar%` else `%do%`
  iter <- iapply(init, 3)
  iterc <- icount(n.runs)
    
  model.out <- foreach(run=iterc, init=iter, .combine=abind2,
                         .multicombine=TRUE) %op% {
    stochastic.ev = eval(stochastic.e)
    SODModelrun(parms, times, locations, init, lambda.ex, K,
                stochastic.e=stochastic.ev, stochastic.d,
                verbose=(verbose & !parallel), run) 
  }

  if (length(dim(model.out))==3) model.out <- abind(model.out, along=4)
  names(dimnames(model.out)) <- c("Time", "Location", "Class", "Replicate")

  if(verbose & !parallel) {
    setTxtProgressBar(z, n.runs*(tail(times, 1) - times[1] + 1))
    close(z)
  }
  
  class(model.out) <- c("SODD", "array")
  return(model.out)
}


#TODO: Make times argument more flexible: scalar or vector, non-consective
#means thinning output
#'Runs the disease model.  Outputs a large matrix of population by species, 
#'sizeclass, location
#'@import plyr reshape2
SODModelrun <- function(parms, times, locations, init, lambda.ex, K,
                        stochastic.e, stochastic.d, verbose, run=NULL) {
  #Tests
  if(!all.equal(dim(init), c(nrow(locations), 2*nrow(parms)))) {
    stop("Dimensions of initial values and parameters do not match")
  }
  
  if(!(length(lambda.ex) == 1 || 
       length(lambda.ex) == length(times) ||
       length(lambda.ex) == nrow(locations) ||
       length(lambda.ex) == nrow(locations) * length(times))) {
    stop("lambda.ex must be 1 or size of match time steps, 
          number of locations, or both")
  }
  
  if(!(is.null(stochastic.e) ||
       (class(stochastic.e) == "numeric" && 
        length(stochastic.e) == length(times)))) {
    stop("stochastic.e must be FALSE, NULL, or a numeric vector of the same
          length as times")
  }
  
  
  n.locations <- nrow(locations)
  
  #Clean up the parms data frame
  
  parms$class <- 1:nrow(parms)
  parms$kernel.fn <- as.character(parms$kernel.fn)
  parms$species <- factor(parms$species, 
                             levels= as.character(unique(parms$species)))
  
  #Convert parameter data frame into list of parameters
  parms.obj <- MakeParmsList(parms)
  
  #Create dispersal matrices
  spread.matrices <- MakeDispMatrix(parms, locations, parms.obj)  
  
  #Unload list of parms into memory
  list2env(parms.obj, envir=environment())

  
  #Create the transition matrix, which includes size transitions, mortality,
  # and resprouting, but not recruitment or infection
  
  tran.mat <- matrix(0, nrow=n.classes*2, ncol=n.classes*2)
  diags <- row(tran.mat) - col(tran.mat)
  fec.mat <- tran.mat
  diag(tran.mat) <- 1 - trans.vec - mort.vec
  tran.mat[diags==2] <- trans.vec[1:(n.classes*2 - 2)]
  
  for(i in 1:n.species) {
    classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
    tran.mat[classindex[1],classindex] <- tran.mat[classindex[1],classindex] + 
                                          (resprout.vec*mort.vec)[classindex] 
  }
  
  # Create empty data matrix and populate with initial values
  
  pop <- array(NA, dim=c(length(times), n.locations, 2*n.classes), 
               dimnames=list(Time=times,
                             Location=1:n.locations, 
                             Class=paste(rep(names(classes),each=2),
                                         rep(c("S","I"),n.classes),sep=".")))
  
  
  
  pop[1,,] <- init
  spore.burden <- matrix(NA, nrow=n.classes+1, ncol=n.locations)

  # Create matrix of external spore burden at all times and locations
  
  if (length(lambda.ex == 1)) spore.burden.ex <- matrix(lambda.ex, nrow=length(times), ncol=n.locations)
  if (length(lambda.ex == length(times))) spore.burden.ex <- matrix(lambda.ex, nrow=length(times), ncol=n.locations)
  if (length(lambda.ex == n.locations)) spore.burden.ex <- matrix(lambda.ex, nrow=length(times), ncol=n.locations, byrow=TRUE)
  if (length(lambda.ex) == length(locations) * length(times)) spore.burden.ex <- lambda.ex
  
  # Create a time vector of environmental stochasticity if none is given
  
  if (is.null(stochastic.e)) {
    stochastic.e <- rep(1,length(times))
  }
    

  
  # Run the simulation
  
  for(time in times[-(length(times))]) {
    
    for(class in classes) {
      spore.burden[class,] <- pop[time,,class*2] %*% spread.matrices[class,,]
    }
    spore.burden[n.classes+1,] <- spore.burden.ex[time,]
    
    force.infection <- (waifw %*% spore.burden) * stochastic.e[time]
    infection.rate <- 1 - exp(-force.infection)
    real.recovery <- (1-infection.rate) * recover
    
    if(stochastic.d==FALSE) {
    
    for(location in 1:n.locations) {
        infect <- infection.rate[,location]
        real.rec <- real.recovery[,location]
        inf.matrix <- diag(as.vector(rbind(1-infect,1-real.rec)))
        inf.matrix[row(inf.matrix) - col(inf.matrix) == 1] <- as.vector(rbind(0,infect))[-1]
        inf.matrix[row(inf.matrix) - col(inf.matrix) == -1] <- as.vector(rbind(0,real.rec))[-1]
        
        pop.inf <- inf.matrix %*% pop[time,location,]
        
        E <- DensityDependence(pop[time,location,], compete, K)
        for(i in 1:n.species) {
          classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
          fec.mat[classindex[1],classindex] <- (pmax(E*recruit.vec,0))[classindex]
        }
        trans.mat <- tran.mat + fec.mat
        pop[time + 1,location,] <- trans.mat %*% pop.inf
      }
    } else {
      
      for(location in 1:n.locations) {
        infect <- infection.rate[,location]
        real.rec <- real.recovery[,location]
        inf.matrix <- diag(as.vector(rbind(1-infect,1-real.rec)))
        inf.matrix[row(inf.matrix) - col(inf.matrix) == 1] <- as.vector(rbind(0,infect))[-1]
        inf.matrix[row(inf.matrix) - col(inf.matrix) == -1] <- as.vector(rbind(0,real.rec))[-1]
        
        pop.inf <- colSums(aaply(1:(n.classes*2),1,function(x) rmultinom(n=1, size=pop[time,location,x], prob=inf.matrix[,x])))
        
        pop.tran <- colSums(aaply(1:(n.classes*2),1, function(x) rmultinom(n=1, size=pop.inf[x], prob=tran.mat[,x])))
        
        E <- DensityDependence(pop[time,location,], compete, K)
        recruitment <- rpois(n=n.classes*2, lambda=pmax(E*pop.inf*recruit.vec,0))  #TODO: Fix this so it doesn't just put a floor, but kills first-year recruits.  Also the pmax above.
        recruits <- rep(0, n.classes*2)
        for(i in 1:n.species) {
          classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
          recruits[classindex[1]] <- sum(recruitment[classindex])
        }
        
        pop[time + 1,location,] <- pop.tran + recruits
      }
    }
  
  if(verbose) setTxtProgressBar(get("z", parent.frame()), time + (length(times-1)*(run-1)))
  }
  

  attr(pop, "spread.matrices") <- spread.matrices
  attr(pop, "locations") <- locations
  return(pop)
  
}

#'Create a list object from the data frame of parameters
#'@import gdata
MakeParmsList <- function(parms) {
  
  classes <- 1:nrow(parms)
  names(classes) <- paste(parms$species, parms$sizeclass, sep=".")
  
  parms$kernel.fn <- as.character(parms$kernel.fn)
  parms$species <- factor(parms$species, 
                             levels= as.character(unique(parms$species)))
  
  sizeclasses <- as.vector(table(parms$species))
  names(sizeclasses) <- 1:length(unique(parms$species))
  
  parms.obj <- list(
    classes = classes,
    sizeclasses = sizeclasses,
    n.species = length(unique(parms$species)),
    n.classes = nrow(parms),
    waifw = cbind(as.matrix(parms[,matchcols(parms, 
                                              with="waifw[0-9]+"),],
                      dimnames = list(names(classes),names(classes))), waifw.ex=parms$waifw.ex),
    recruit.vec = as.vector(rbind(parms$S.recruit, 
                                  parms$I.recruit)),
    mort.vec = as.vector(rbind(parms$S.mortality, 
                               parms$I.mortality)),
    trans.vec = as.vector(rbind(parms$S.transition, 
                                parms$I.transition)),
    resprout.vec = as.vector(rbind(parms$S.resprout, 
                                   parms$I.resprout)),
    recover = parms$recover,
    compete = parms$compete
  )
  return(parms.obj)
}

#'Generate a dispersal matrix from locations.  The matrix calculates pairwise
#'distances between locations weighted by the dispersal kernel of each class
#'@import plyr
MakeDispMatrix <- function(parms, locations, parms.obj) {
  
  distance.matrix <- as.matrix(dist(locations))
  
  # Generate dispersal matrices for each class
  spread.matrices <- daply(parms,.(class), function(x) {
    do.call(x$kernel.fn, 
            unname(c(list(distance.matrix), 
                     as.list(na.omit(x[matchcols(x, "kernel.par[0-9]+")]))
            ))
    )
  })
  return(spread.matrices)
}  


# First plot: Big and small tanoaks, mean and many replicates on the same plot

