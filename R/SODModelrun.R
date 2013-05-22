#TODO: Make times argument more flexible: scalar or vector, non-consective
#means thinning output
#'Runs the disease model.  Outputs a large matrix of population by species, 
#'sizeclass, location
#'@import plyr reshape2
SODModelrun <- function(parms, times, locs, init, spread.matrices, lambda.ex,
                        stochastic.e, stochastic.d, verbose, run=NULL) {
  list2env(parms, envir=environment())

#   #Tests
#   if(!all.equal(dim(init), c(nrow(locations), 2*nrow(parms)))) {
#     stop("Dimensions of initial values and parameters do not match")
#   }
#   
#   if(!(length(lambda.ex) == 1 || 
#        length(lambda.ex) == length(times) ||
#        length(lambda.ex) == nrow(locations) ||
#        length(lambda.ex) == nrow(locations) * length(times))) {
#     stop("lambda.ex must be 1 or size of match time steps, 
#           number of locations, or both")
#   }
#   
#   if(!(is.null(stochastic.e) ||
#        (class(stochastic.e) == "numeric" && 
#         length(stochastic.e) == length(times)))) {
#     stop("stochastic.e must be FALSE, NULL, or a numeric vector of the same
#           length as times")
#   }
  
  
  
  #Clean up the parms data frame
  
  #Convert parameter data frame into list of parameters
  
  #Unload list of parms into memory

  
  #Create the transition matrix, which includes size transitions, mortality,
  # and resprouting, but not recruitment or infection
  
  tran.mat <- matrix(0, nrow=n.class*2, ncol=n.class*2)
  diags <- row(tran.mat) - col(tran.mat)
  fec.mat <- tran.mat
  diag(tran.mat) <- 1 - trans.vec - mort.vec
  tran.mat[diags==2] <- trans.vec[1:(n.class*2 - 2)]
  
  for(i in 1:n.species) {
    classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
    tran.mat[classindex[1],classindex] <- tran.mat[classindex[1],classindex] + 
                                          (resprout.vec*mort.vec)[classindex] 
  }
  
  # Create empty data matrix and populate with initial values
  
  pop <- array(NA, dim=c(length(times), n.loc, 2*n.class)) 
  
  pop[1,,] <- init
  spore.burden <- matrix(NA, nrow=n.class+1, ncol=n.loc)

  # Create matrix of external spore burden at all times and locations
  
  if (length(lambda.ex == 1)) spore.burden.ex <- matrix(lambda.ex, nrow=length(times), ncol=n.loc)
  if (length(lambda.ex == length(times))) spore.burden.ex <- matrix(lambda.ex, nrow=length(times), ncol=n.loc)
  if (length(lambda.ex == n.loc)) spore.burden.ex <- matrix(lambda.ex, nrow=length(times), ncol=n.loc, byrow=TRUE)
  if (length(lambda.ex) == length(locs) * length(times)) spore.burden.ex <- lambda.ex
  
  # Create a time vector of environmental stochasticity if none is given
  
  if (is.null(stochastic.e)) {
    stochastic.e <- rep(1,length(times))
  }
    

  
  # Run the simulation
  
  for(time in times[-(length(times))]) {
    
    for(class in 1:n.class) {
      spore.burden[class,] <- pop[time,,class*2] %*% spread.matrices[class,,]
    }
    spore.burden[n.class+1,] <- spore.burden.ex[time,]
    
    force.infection <- (waifw %*% spore.burden) * stochastic.e[time]
    infection.rate <- 1 - exp(-force.infection)
    real.recovery <- (1-infection.rate) * recover
    
    if(stochastic.d==FALSE) {
    
    for(loc in 1:n.loc) {
        infect <- infection.rate[,loc]
        real.rec <- real.recovery[,loc]
        inf.matrix <- diag(as.vector(rbind(1-infect,1-real.rec)))
        inf.matrix[row(inf.matrix) - col(inf.matrix) == 1] <- as.vector(rbind(0,infect))[-1]
        inf.matrix[row(inf.matrix) - col(inf.matrix) == -1] <- as.vector(rbind(0,real.rec))[-1]
        
        pop.inf <- inf.matrix %*% pop[time,loc,]
        
        E <- DensityDependence(pop[time,loc,], compete, K)
        for(i in 1:n.species) {
          classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
          fec.mat[classindex[1],classindex] <- (pmax.int(E*recruit.vec,0))[classindex]
        }
        trans.mat <- tran.mat + fec.mat
        pop[time + 1,loc,] <- trans.mat %*% pop.inf
      }
    } else {
      
      for(loc in 1:n.loc) {
        infect <- infection.rate[,loc]
        real.rec <- real.recovery[,loc]
        inf.matrix <- diag(as.vector(rbind(1-infect,1-real.rec)))
        inf.matrix[row(inf.matrix) - col(inf.matrix) == 1] <- as.vector(rbind(0,infect))[-1]
        inf.matrix[row(inf.matrix) - col(inf.matrix) == -1] <- as.vector(rbind(0,real.rec))[-1]
        
        pop.inf <- colSums(aaply(1:(n.class*2),1,function(x) rmultinom(n=1, size=pop[time,loc,x], prob=inf.matrix[,x])))
        
        pop.tran <- colSums(aaply(1:(n.class*2),1, function(x) rmultinom(n=1, size=pop.inf[x], prob=tran.mat[,x])))
        
        E <- DensityDependence(pop[time,loc,], compete, K)
        recruitment <- rpois(n=n.class*2, lambda=pmax.int(E*pop.inf*recruit.vec,0))  #TODO: Fix this so it doesn't just put a floor, but kills first-year recruits.  Also the pmax above.
        recruits <- rep(0, n.class*2)
        for(i in 1:n.species) {
          classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
          recruits[classindex[1]] <- sum(recruitment[classindex])
        }
        
        pop[time + 1,loc,] <- pop.tran + recruits
      }
    }
  
  if(verbose) setTxtProgressBar(get("z", parent.frame()), time + (length(times-1)*(run-1)))
  }
  return(pop)
  
}

#'Create a list object from the data frame of parameters
MakeParmsList <- function(parms, locs) {
  parms.obj <- within(as.list(parms), {
    n.loc <- nrow(locs)
    loc.names <- rownames(locs)
    coord.names <- colnames(locs)
    n.species = length(unique(parms$species))
    n.class <- nrow(parms)
    sizeclasses <- as.vector(table(species))
    class.names <- paste0(species, ".", sizeclass)
    K <- ifelse(exists("K"), K[1], 1)

    waifw <- do.call(cbind, c(mget(ls(pattern="waifw[0-9]+")),
                              list(waifw.ex=waifw.ex)))/K
    species <- factor(species, levels= as.character(unique(parms$species)))
    kernel.fn <- as.character(kernel.fn)
    kernel.pars <- alply(do.call(cbind, mget(ls(pattern="kernel.par[0-9]"))),
                         1, na.omit)
    recruit.vec = as.vector(rbind(S.recruit, I.recruit))
    mort.vec <- as.vector(rbind(S.mortality, I.mortality))
    trans.vec <- as.vector(rbind(S.transition,I.transition))
    resprout.vec <- as.vector(rbind(S.resprout,I.resprout))
    rm(list=c("waifw.ex", ls(pattern="waifw[0-9]+"),
              ls(pattern="kernel.par[0-9]")))
  })
  parms.obj <- parms.obj[!sapply(parms.obj, is.null)]
  return(parms.obj)
}

#'Generate a dispersal matrix from locations.  The matrix calculates pairwise
#'distances between locations weighted by the dispersal kernel of each class
#'@import plyr
MakeDispMatrix <- function(parms, locs) {
  
  distance.matrix <- as.matrix(dist(locs))
  
  # Generate dispersal matrices for each class
  spread.matrices <- aaply(1:n.class, 1, function(x) {
    do.call(parms$kernel.fn[x],
            unname(c(list(distance.matrix), kernel.pars[[x]])))
  }, .drop=FALSE)
  dimnames(spread.matrices) <- list(Class=parms$class.names,
                                    parms$loc.names, parms$loc.names)
  return(spread.matrices)
}  


# First plot: Big and small tanoaks, mean and many replicates on the same plot
