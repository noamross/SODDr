#'Runs the disease model.  Outputs a large matrix of population by species, ageclass, location
#'@import plyr 
#'@importFrom tidyr separate_
#'@export
SODModel <- function(parms.df, locations, time.steps, init, df.out=TRUE, verbose=interactive()) {
  
  #Tests
  if(!all.equal(dim(init), c(nrow(locations), 2*nrow(parms.df)))) {
    stop("Dimensions of initial values and parameters do not match")
  }
  
  parms.df$class <- 1:nrow(parms.df)
  n.locations <- nrow(locations)
  
  #Convert parameter data frame into list of parameters
  parms.obj <- MakeParmsList(parms.df)
  #Create dispersal matrices
  spread.matrices <- MakeDispMatrix(parms.df, locations, parms.obj)  
  
  # Create transition matrix
  # TODO: Make this into a function
  
  #Unload list of parms into memory
  for(i in 1:length(parms.obj)) assign(names(parms.obj)[i], parms.obj[[i]])
  
  tran.mat <- matrix(0, nrow=n.classes*2, ncol=n.classes*2)
  diags <- row(tran.mat) - col(tran.mat)
  fec.mat <- tran.mat
  diag(tran.mat) <- 1 - trans.vec - mort.vec
  
  tran.mat[diags==1] <- diag(tran.mat)[-(nrow(tran.mat))] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
  tran.mat[diags==-1] <- diag(tran.mat)[-1] * rep(c(1,0),n.classes)[-(nrow(tran.mat))]
  
  tran.mat[diags==2] <- trans.vec[1:(n.classes*2 - 2)]
  tran.mat[diags==3] <- trans.vec[1:(n.classes*2 -3)] * rep(c(1,0),n.classes)[1:(n.classes*2 - 3)]
  tran.mat[diags==1] <- tran.mat[diags==1] + trans.vec[1:(n.classes*2-1)] * rep(c(0,1), n.classes)[1:(n.classes*2 - 1)]
  
  # Step update is (trans.mat * force.mat + fec.mat) * population
  
  # Create empty data matrix and populate with initial values
  
  pop <- array(NA, dim=c(length(time.steps), n.locations, 2*n.classes), 
               dimnames=list(Time=time.steps,
                             Location=1:n.locations, 
                             Class=paste(rep(names(classes),each=2),
                                         rep(c("S","I"),n.classes),sep=",")))
  
  
  
  pop[1,,] <- init
  spore.burden <- matrix(NA, nrow=n.classes, ncol=n.locations)
  #Rprof("out.prof")

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
    
    force.infection <- waifw %*% spore.burden
    real.recovery <- (1-force.infection) * recover
    
    for(location in 1:n.locations) {
      force <- force.infection[,location]
      real.rec <- real.recovery[,location]
      force.matrix <- matrix(rbind(c(1-force, force), c(real.rec, 1-real.rec)), 2*n.classes, 2*n.classes, byrow=TRUE)
      
      E <- DensityDependence(pop[time,location,], space)
      for(i in 1:n.species) {
        classindex <- 1:(2*ageclasses[i]) + (sum(ageclasses[0:(i-1)])*2)
        fec.mat[classindex[1],classindex] <- (E*recruit.vec + resprout.vec*mort.vec)[classindex]  #Fecundities + Death X Resprout probabilties
      }
      trans.mat <- tran.mat*force.matrix + fec.mat
      pop[time + 1,location,] <- trans.mat %*% pop[time,location,]
    }
  if(verbose) z$step()
  }
  if(df.out) {
    pop <- melt(pop, value.name="Population")
    pop <- tidyr::separate_(pop, "Class", into=c("Species", "AgeClass", "Disease"), sep=",", remove=TRUE, convert=TRUE)
  }
  return(pop)
  
}


#'@import gdata
MakeParmsList <- function(treeparms.df) {
  
  classes <- 1:nrow(treeparms.df)
  names(classes) <- paste(treeparms.df$species, treeparms.df$ageclass, sep=",")
  
  ageclasses <- as.vector(table(treeparms.df$species))
  names(ageclasses) <- 1:length(unique(treeparms.df$species))
  
  parms.obj <- list(
    classes = classes,
    ageclasses = ageclasses,
    n.species = length(unique(treeparms.df$species)),
    n.classes = nrow(treeparms.df),
    waifw = as.matrix(treeparms.df[,matchcols(treeparms.df, 
                                              with="waifw[0-9]+"),],
                      dimnames = list(names(classes),names(classes))),
    recruit.vec = as.vector(rbind(treeparms.df$S.recruit, 
                                  treeparms.df$I.recruit)),
    mort.vec = as.vector(rbind(treeparms.df$S.mortality, 
                               treeparms.df$I.mortality)),
    trans.vec = as.vector(rbind(treeparms.df$S.transition, 
                                treeparms.df$I.transition)),
    resprout.vec = as.vector(rbind(treeparms.df$S.resprout, 
                                   treeparms.df$I.resprout)),
    recover = treeparms.df$recover,
    space = treeparms.df$space
  )
  return(parms.obj)
}




#'Generate a dispersal matrix from locations.
#'@import plyr
MakeDispMatrix <- function(parms.df, locations, parms.obj) {
  
  distance.matrix <- as.matrix(dist(locations[,c("x","y")]))
  
  # Generate dispersal matrices for each class
  spread.matrices <- daply(parms.df,"class", function(x) {
    do.call(x$kernel.fn, 
            unname(c(list(distance.matrix), 
                     as.list(na.omit(x[matchcols(x, "kernel.par[0-9]+")]))
            ))
    )
  })
  return(spread.matrices)
}
  