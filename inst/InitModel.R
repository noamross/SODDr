#' Generate a Model from which to draw a random forest based on plot data
#' @import sp mgcv bbmle automap gstat
#' @export
InitModel <- function(data, plot.area=500, 
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
  #Fit a model to the data
  if(!is.null(form) & length(model.type)==1) {
    plot.model <- do.call(model.type, list(formula=form, data=plot.data,
                       family=poisson()))
  } else {
    models <- list()
    for(mt in model.type) {
      if(mt=="gam") {
        frms <- c(as.formula(paste0(Count,"~1")),
                  as.formula(paste0(Count, "~s(", xy[1], ")")),
                  as.formula(paste0(Count, "~s(", xy[2], ")")),
                  as.formula(paste0(Count, "~s(", xy[1], ") + s(", xy[2],")")),
                  as.formula(paste0(Count, "~s(", xy[1], ",", xy[2], ")")))
      } else {
        frms <- c(as.formula(paste0(Count,"~1")),
                  as.formula(paste0(Count, "~", xy[1])),
                  as.formula(paste0(Count, "~", xy[2])),
                  as.formula(paste0(Count, "~", xy[1], "+", xy[2])),
                  as.formula(paste0(Count, "~", xy[1], "*", xy[2])))
      }
      for(i in 1:length(frms)) {
        models[[paste0(mt,i)]] <- try(
                                      do.call(mt, list(formula=frms[[i]],
                                              data=as.data.frame(plot.data),
                                              family=poisson())),
                                      silent=TRUE)
      }
    }
    models <- models[which(!(llply(models,class)=="try-error"))]
    model.order <- rank(laply(models,AIC))
    for (m in model.order) {
      plot.model <- models[[m]]
      model.plot.predictions <- try(predict(plot.model, type="link"))
      model.grid.predictions <- try(predict(plot.model, newdata=map.model, 
                                            type="link", se.fit=TRUE))
      if(!(class(model.plot.predictions) == "try-error" | 
             class(model.grid.predictions) == "try-error")) break
    }
  }

  #Generate predictions and residuals
  plot.data$residuals <- residuals(plot.model, type="working")
  #Generate a kriging surface on the residuals and apply to the hex grid
  #krig <- autoKrige(residuals~1, plot.data, new_data=hex.pts)$krige_output
  krig <- gstat::idw(residuals~1, plot.data, newdata=hex.pts, debug.level=0)
  attributes(model.grid.predictions$fit) <- NULL
  attributes(model.grid.predictions$se.fit) <- NULL
  map.model$mod.pred <- model.grid.predictions$fit
  map.model$krig.fit <- krig$var1.pred
  map.model$mean <- krig$var1.pred + model.grid.predictions$fit
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
  
  
  