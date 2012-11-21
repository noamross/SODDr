#' Spore dispersal via normal distribution
#' Assumes a regular grid
#' 
#' By default, the normal distribution is 
#' @param var
#' @param gridsize
#' @param varcov
require(cubature)
require(mnormt)

adaptIntegrate(
dmnorm(
  
normal.dispersal <- function(distance, var=1, gridsize=1, varcov=NULL)
  
  if(is.null(varcov)) varcov <- var*rbind(c(var, 0), c(0, var))

  if(distance==0)
    
Step 1: Run model with point-based mvtnorm kernel (should include option for normalizing to zero)
Step 2: Step two - divide dispersal options into point-based or area-based (integration)
Step 3: 
  
if(!(class(locations)=="lattice") || is.numeric(attr(locations, "gridsize"))) {
  stop("Integrated dispersal methods require regular grid location matrices")
}

gridsize <- attr(locations, "gridsize")

spread.matrix <- matrix(NA, nrow(locations), nrow(locations))
spread.matrix2 <- spread.matrix
for(i in 1:nrow(locations)) {
    for(j in 1:nrow(locations)) {
      spread.matrix[i,j] <- sadmvn(lower=c(locations[j, "x"] - 0.5*gridsize,
                                           locations[j, "y"] - 0.5*gridsize),
                                   upper=c(locations[j, "x"] + 0.5*gridsize,
                                           locations[j, "y"] + 0.5*gridsize),
                                   mean=c(locations[i, "x"], locations[i, "y"]),
                                   varcov=varcov)
    }
}
spread.matrix2 <- dnorm(as.matrix(dist(locations[,2:3])))
    

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