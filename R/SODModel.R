#'Runs the disease model.  Outputs a large matrix of population by species, sizeclass, location
#'@import plyr reshape2
#'@export
SODModel <- function(parms.df, locations, time.steps, init, lambda.ex = 0, df.out=TRUE, verbose=interactive(), stochastic.d=FALSE, stochastic.e=FALSE) {
  
  #Tests
  if(!all.equal(dim(init), c(nrow(locations), 2*nrow(parms.df)))) {
    stop("Dimensions of initial values and parameters do not match")
  }
  
  if(!(length(lambda.ex) == 1 || 
       length(lambda.ex) == length(time.steps) ||
       length(lambda.ex) == length(locations) ||
       length(lambda.ex) == length(locations) * length(time.steps))) {
    stop("lambda.ex must be 1 or size to match time steps, 
          number of locations, or both")
  }
  
  if(!(stochastic.e == FALSE ||
       stochastic.e == NULL ||
       (class(stochastic.e) == "numeric" & 
        length(stochastic.e) == length(time.steps)))) {
    stop("stochastic.e must be FALSE, NULL, or a numeric vector of the same
          length as time.steps")
  }
  
  
  n.locations <- nrow(locations)
  
  #Clean up the parms data frame
  
  parms.df$class <- 1:nrow(parms.df)
  parms.df$kernel.fn <- as.character(parms.df$kernel.fn)
  parms.df$species <- factor(parms.df$species, 
                             levels= as.character(unique(parms.df$species)))
  
  #Convert parameter data frame into list of parameters
  parms.obj <- MakeParmsList(parms.df)
  #Create dispersal matrices
  spread.matrices <- MakeDispMatrix(parms.df, locations, parms.obj)  
  
  # Create transition matrix
  # TODO: Make this into a function
  
  #Unload list of parms into memory
  list2env(parms.obj, envir=environment())
#  for(i in 1:length(parms.obj)) assign(names(parms.obj)[i], parms.obj[[i]])
  
  tran.mat <- matrix(0, nrow=n.classes*2, ncol=n.classes*2)
  diags <- row(tran.mat) - col(tran.mat)
  fec.mat <- tran.mat
  diag(tran.mat) <- 1 - trans.vec - mort.vec
  
#  tran.mat[diags==1] <- diag(tran.mat)[-(nrow(tran.mat))] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
#  tran.mat[diags==-1] <- diag(tran.mat)[-1] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
  
  tran.mat[diags==2] <- trans.vec[1:(n.classes*2 - 2)]
#  tran.mat[diags==3] <- trans.vec[1:(n.classes*2 -3)] * rep(c(1,0),n.classes)[1:(n.classes*2 - 3)]
#  tran.mat[diags==1] <- tran.mat[diags==1] + trans.vec[1:(n.classes*2-1)] * rep(c(0,1), n.classes)[1:(n.classes*2 - 1)]
  for(i in 1:n.species) {
    classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
    tran.mat[classindex[1],classindex] <- tran.mat[classindex[1],classindex] + (resprout.vec*mort.vec)[classindex]  #Fecundities + Death X Resprout probabilties
  }
  # Step update is (trans.mat * force.mat + fec.mat) * population
  
  # Create empty data matrix and populate with initial values
  
  pop <- array(NA, dim=c(length(time.steps), n.locations, 2*n.classes), 
               dimnames=list(Time=time.steps,
                             Location=1:n.locations, 
                             Class=paste(rep(names(classes),each=2),
                                         rep(c("S","I"),n.classes),sep=",")))
  
  
  
  pop[1,,] <- init
  spore.burden <- matrix(NA, nrow=n.classes+1, ncol=n.locations)
  #Rprof("out.prof")

  #Create matrix of external spore burden
  
  if (length(lambda.ex == 1)) spore.burden.ex <- matrix(lambda.ex, nrow=length(time.steps), ncol=n.locations)
  if (length(lambda.ex == length(time.steps))) spore.burden.ex <- matrix(lambda.ex, nrow=length(time.steps), ncol=n.locations)
  if (length(lambda.ex == n.locations)) spore.burden.ex <- matrix(lambda.ex, nrow=length(time.steps), ncol=n.locations, byrow=TRUE)
  if (length(lambda.ex) == length(locations) * length(time.steps)) spore.burden.ex <- lambda.ex
  
  if(verbose) {
  oldopt <- getOption("warn")
  options(warn=2)
  z <- try(create_progress_bar("time"), silent=TRUE)
  if(class(z)=="try-error") z <- try(create_progress_bar("text"))
  options(warn=oldopt)
  z$init(length(time.steps)-1)
  on.exit(z$term)
  }
  for(time in time.steps[-(length(time.steps))]) {
    
    # First act in simulation step.  Given population at each location, calculate spore burden at each location
    for(class in classes) {
      spore.burden[class,] <- pop[time,,class*2] %*% spread.matrices[class,,]
    }
    spore.burden[n.classes+1,] <- spore.burden.ex[time,]
    
    force.infection <- waifw %*% spore.burden
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
        
        E <- DensityDependence(pop[time,location,], compete)
        for(i in 1:n.species) {
          classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
          fec.mat[classindex[1],classindex] <- (E*recruit.vec)[classindex]
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
        
        E <- DensityDependence(pop[time,location,], compete)
        recruitment <- rpois(n=n.classes*2, lambda=E*pop.inf*recruit.vec)
        recruits <- rep(0, n.classes*2)
        for(i in 1:n.species) {
          classindex <- 1:(2*sizeclasses[i]) + (sum(sizeclasses[0:(i-1)])*2)
          recruits[classindex[1]] <- sum(recruitment[classindex])
        }
        
        pop[time + 1,location,] <- pop.tran + recruits
      }
    }
  
  if(verbose) z$step()
  }
  if(df.out) {
    pop.df <- melt(pop, value.name="Population")
    pop.df$Class <- factor(pop.df$Class, levels <- dimnames(pop)$Class)
    pop.df <- arrange(pop.df, Time,Location,Class)
    pop.df$Species <- factor(rep(parms.df$species, each=2)) #, levels=levels(parms.df$species))
    pop.df$SizeClass <- factor(rep(parms.df$sizeclass, each=2))
    pop.df$Disease <- factor(c("S","I"), c("S","I"))
    pop.df$Class <- NULL
    pop.df <- pop.df[,c(1,2,4,5,6,3)]
    attr(pop.df, "spread.matrices") <- spread.matrices
    return(pop.df)
  }
  attr(pop, "spread.matrices") <- spread.matrices
  return(pop)
  
}


#'@import gdata
MakeParmsList <- function(parms.df) {
  
  classes <- 1:nrow(parms.df)
  names(classes) <- paste(parms.df$species, parms.df$sizeclass, sep=".")
  
  parms.df$kernel.fn <- as.character(parms.df$kernel.fn)
  parms.df$species <- factor(parms.df$species, 
                             levels= as.character(unique(parms.df$species)))
  
  sizeclasses <- as.vector(table(parms.df$species))
  names(sizeclasses) <- 1:length(unique(parms.df$species))
  
  parms.obj <- list(
    classes = classes,
    sizeclasses = sizeclasses,
    n.species = length(unique(parms.df$species)),
    n.classes = nrow(parms.df),
    waifw = cbind(as.matrix(parms.df[,matchcols(parms.df, 
                                              with="waifw[0-9]+"),],
                      dimnames = list(names(classes),names(classes))), waifw.ex=parms.df$waifw.ex),
    recruit.vec = as.vector(rbind(parms.df$S.recruit, 
                                  parms.df$I.recruit)),
    mort.vec = as.vector(rbind(parms.df$S.mortality, 
                               parms.df$I.mortality)),
    trans.vec = as.vector(rbind(parms.df$S.transition, 
                                parms.df$I.transition)),
    resprout.vec = as.vector(rbind(parms.df$S.resprout, 
                                   parms.df$I.resprout)),
    recover = parms.df$recover,
    compete = parms.df$compete
  )
  return(parms.obj)
}




#'Generate a dispersal matrix from locations.  The matrix calculates pairwise
#'distances between locations weighted by the dispersal kernel of each class
#'@import plyr
MakeDispMatrix <- function(parms.df, locations, parms.obj) {
  
  distance.matrix <- as.matrix(dist(locations[,c("x","y")]))
  
  # Generate dispersal matrices for each class
  spread.matrices <- daply(parms.df,.(class), function(x) {
    do.call(x$kernel.fn, 
            unname(c(list(distance.matrix), 
                     as.list(na.omit(x[matchcols(x, "kernel.par[0-9]+")]))
            ))
    )
  })
  return(spread.matrices)
}
  