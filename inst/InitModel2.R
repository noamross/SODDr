#' Generate a Model from which to draw a random forest based on plot data
#' @import sp mgcv bbmle automap gstat
#' @export
InitModel2 <- function(data, plot.area=500, 
                      model.type=c("glm","gam"), Count="Count",
                      form=NULL, xy = c("Easting", "Northing"),
                      p4string="+proj=utm +zone=10 +ellps=intl +units=m") {

  # Make plot.data into a SpatialPointsDataFrame
  if(class(data) == "data.frame") {
    coordinates(data) <- xy
    proj4string(data) <- CRS(p4string)
  }

  #Set up boundaries of the hexagonal grid
  bb <- bbox(data)
  
  #Make a SpatialPoints hexagonal grid in the bounding box
  hex.pts <- spsample(data, type = "hexagonal", 
                     cellsize = sqrt(plot.area*2/(sqrt(3))), bb=bb)
  locations <- coordinates(hex.pts)
  rownames(locations) <- 1:nrow(locations)
  colnames(locations) <- xy
  map.model <- as.data.frame(hex.pts)
  names(map.model) <- xy
  fnenvt <- environment()
  formulas <- list()
  modelist <- list()
  o <- ddply(as.data.frame(data), c("Species","SizeClass", "Infected"), function(plot.data) {
  
  
  coordinates(plot.data) <- xy
  proj4string(plot.data) <- CRS(p4string)
  plot.data$x <- coordinates(plot.data)[,1]
  plot.data$y <- coordinates(plot.data)[,2]
  #Fit a model to the data
  ak <- autoKrige(Count~x*y, input_data=plot.data, new_data=hex.pts)
  
  coords <- coordinates(plot.data)
  g1 <- glm(Count~1, family=poisson, data=plot.data)
  beta.starting <- coef(g1)
  beta.tuning <- t(chol(vcov(g1)))
  
  n.batch <- 1000
  batch.length <- 100
  n.samples <- n.batch*batch.length
  
  m.1 <- spGLM(z~x+y, family="poisson", coords=coords,
             starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             priors=list("beta.Flat", "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             tuning=list("beta"=c(0.1,0.1,0.1), "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)
  
  burn.in <- 0.9*n.samples
  sub.samps <- burn.in:n.samples
  
  plot(m.1$p.beta.theta.samples)
  
  beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]
  x <- as.matrix(rep(1,n))
y.hat <- apply(exp(x%*%beta.hat+w.hat), 2, function(x){rpois(n, x)})


  
  require(MBS)
  
  #Generate predictions and residuals
  #plot.data$residuals <- residuals(plot.model, type="working")
  #Generate a kriging surface on the residuals and apply to the hex grid
  #krig <- autoKrige(residuals~1, plot.data, new_data=hex.pts)$krige_output
  #krig <- gstat::idw(residuals~1, plot.data, newdata=hex.pts, debug.level=0)
  #attributes(model.grid.predictions$fit) <- NULL
  #attributes(model.grid.predictions$se.fit) <- NULL
  #map.model$mod.pred <- model.grid.predictions$fit
  #map.model$krig.fit <- krig$var1.pred
  map.model$mean <- ak$var1.pred
  map.model$sd <- model.grid.predictions$se.fit
  map.model$expected <- exp(map.model$mean)
  map.model$location <- 1:nrow(map.model)
  model.formulas <- formula
  assign("formulas", append(get("formulas", envir=fnenvt), formula(plot.model)), envir=fnenvt)
  assign("modelist", c(get("modelist", envir=fnenvt), list(summary(plot.model))), envir=fnenvt)
  return(map.model)
  })
  o <- SpatialPointsDataFrame(coords=o[,c(4,5)], data=o,coords.nrs=c(4,5),proj4string=CRS(p4string))
  attr(o, "formulas") <- formulas
  attr(o, "modelist") <- modelist
  attr(o, "locations") <- locations
  return(o)
}

#' @import reshape2
#' @export
InitDraw <- function(map.model, error=TRUE, array=TRUE) {
  if(error) {
    map.model$SimCount <- rpois(n=nrow(map.model), 
                                lambda=exp(rnorm(nrow(map.model), 
                                                 map.model$mean, 
                                                 sd=map.model$sd)))
  } else {
    map.model$SimCount <- rpois(n=nrow(map.model), lambda=map.model$expected)
  }                                                 
  if(array) {
    map.model <- as.data.frame(map.model)
    keep.names <- unique(paste(map.model$Species, map.model$SizeClass, 
                               map.model$Infected, sep="_"))
    map.model <- dcast(map.model, location ~ Species + SizeClass + Infected, 
                       value.var="SimCount", drop=FALSE)
    map.model <- map.model[which(names(map.model) %in% keep.names)]
    map.model <- as.matrix(map.model)

  }
  return(map.model)
}
  
  
  