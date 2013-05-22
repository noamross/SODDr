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
SODModel <- function(parms, times, locs, init, reps=1, lambda.ex = 0,
                     stochastic.e=NULL, stochastic.d=FALSE, 
                     verbose=interactive(), parallel=FALSE) {

  # Run tests on inputs
  #test_SOD_inputs()
  # Generate metadata
  #metadata <- make_metadata(parms, times, loc, init)
  if(class(parms)=="data.frame") parms <- MakeParmsList(parms, locs)
  message("Creating spread matrices...")
  spread.matrices <- MakeDispMatrix(parms, locs)  

  
  init <- abind(rep(list(init),each=reps), along=3)
  n.runs <- dim(init)[3]


  
  `%op%` <- if (parallel) `%dopar%` else `%do%`
  iter <- iapply(init, 3)
  iterc <- icount(n.runs)
    
  message("Initializing loop...")
  if(verbose) {
    message(n.runs, " Simulations of ",
            tail(times, 1) - times[1] + 1, " time steps each")
    if(!parallel) {
      z <- plyr:::txtTimerBar((tail(times, 1) - times[1] + 1)*n.runs)
    }
  }
  model.out <- foreach(run=iterc, init=iter, .combine=abind2,
                         .multicombine=TRUE) %op% {
    stochastic.ev = eval(stochastic.e)
    SODModelrun(parms, times, locs, init, spread.matrices, lambda.ex,
                stochastic.e=stochastic.ev, stochastic.d,
                verbose=(verbose & !parallel), run) 
  }

  if (length(dim(model.out))==3) model.out <- abind(model.out, along=4)
  dimnames(model.out) <- list(Time=times, Locations=parms$loc.names,
                               Class=paste(rep(parms$class.names,each=2),
                                           rep(c("S","I"),parms$n.class),
                                           sep="."),
                               Replicate=1:(dim(model.out)[4]))

  if(verbose & !parallel) {
    setTxtProgressBar(z, n.runs*(tail(times, 1) - times[1] + 1))
    close(z)
  }
  
  class(model.out) <- c("SODD", "array")
  attr(model.out, "spread.matrices") <- spread.matrices
  attr(model.out, "locations") <- locs
  attr(model.out, "parms") <- parms
  return(model.out)
}


