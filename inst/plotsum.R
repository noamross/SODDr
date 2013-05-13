require(plyr)
require(ggplot2)
require(sp)
require(hexbin)
require(grid)
require(SODDr)
require(compiler)
require(abind)
require(spBayes)
require(iterators)
require(foreach)
require(doRNG)
require(doParallel)
registerDoParallel()
require(abind)
# enableJIT(3)

# Set some essential variables
SUBSITE <- "JLSP"  # Which subsite of the plot data 
AREA <- 500

# Load the data and merge subsites names
rawdata <- read.csv("data/redwood_plot_network_full.csv", comment.char="#")
subsites <- read.csv("data/plot-site list.csv")
names(subsites)[2] <- "Subsite"
rawdata <- merge(rawdata, subsites)
rawdata$Infected <- factor(ifelse(rawdata$Infected,"I","S"), levels=c("S","I"))

### Data Conversion ###

# Get initial census of live trees at the chosen subsite
census.2002 <- subset(rawdata, Year==2002 & Mortality==0 & Subsite==SUBSITE)

# Aggregate non-sporulating species and assume no infection in these
sporulators <- c("UMCA", "LIDE")
census.2002$Species <- factor(census.2002$Species, 
                              levels = c(levels(census.2002$Species), "OTH"))
census.2002$Species[which(!(census.2002$Species %in% sporulators))] <- "OTH"
census.2002$Infected[which(!(census.2002$Species %in% sporulators))] <- "S"
census.2002 <- droplevels(census.2002)

# Assign size classes to only the LIDE trees, giving the rest SizeClass=1

LIDE.breaks <- c(0,2,10,30, Inf)
census.2002$SizeClass <- factor(1, levels=as.character(1:4))
census.2002$SizeClass[which(census.2002$Species=="LIDE")] <- 
  cut(census.2002$DBH[which(census.2002$Species=="LIDE")], 
      breaks=LIDE.breaks, labels=as.character(1:4))

# Summarize the data by plot
plot.sum <- ddply(census.2002, c("Plot", "Species", "SizeClass", "Infected"), 
                  summarize, Count=length(Species), .drop=FALSE)

# Eliminate unneccessary classes
plot.sum <- subset(plot.sum, !(Species !="LIDE" & SizeClass !="1"))

#Add the Site data back in
plot.sum <- merge(unique(census.2002[,c("Site", "Subsite", "Plot",
                                        "Easting", "Northing")]), plot.sum)


plot.glms <- dlply(plot.sum, .(Species, SizeClass, Infected), function(x) {
                     glm(Count ~ 1, family="poisson", data=x)
                    })


n.batch <- 100
batch.length <- 50
n.samples <- n.batch*batch.length
n.inits <- 5
burn.in <- ceiling(0.9*n.samples)
sub.samps <- (burn.in+1):n.samples

spGLMs <- llply(plot.glms, function(GLM) {
                  spGLM(Count~1, family="poisson", data=GLM$data,
                  coords=as.matrix(GLM$data[,c("Easting", "Northing")]),
        starting=list("beta"=coef(GLM),"phi"=0.06,"sigma.sq"=vcov(GLM), "w"=0),
        priors=list("beta.Flat", "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
        tuning=list("beta"=0.1, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
        amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
        cov.model="exponential", verbose=FALSE)
                })

plot.sum.sp <- plot.sum
coordinates(plot.sum.sp) <- c("Easting","Northing")
proj4string(plot.sum.sp) <- "+proj=utm +zone=10 +ellps=intl +units=m"
bb <- bbox(plot.sum.sp)
hex.pts <- coordinates(spsample(plot.sum.sp, type = "hexagonal", 
                    cellsize = sqrt(AREA*2/(sqrt(3))), bb=bb))
colnames(hex.pts) <- c("Easting", "Northing")
sPredicts <- llply(spGLMs, function(SPG) {
                     spPredict(SPG, pred.coords=hex.pts, 
                               pred.covars=matrix(1, nrow(hex.pts)),
                               start=burn.in, end=n.samples,
                               thin=(n.samples-burn.in)/n.inits,
                               verbose=FALSE)
                  })

inits <- laply(sPredicts, function(SPR) {
                apply(SPR$p.y.predictive.sample,1, function(x) {
                        rpois(length(x), x) })
                })
inits <- aperm(inits, c(3,1,2))
dimnames(inits) <- list(Location=1:(dim(inits)[1]), Class=names(plot.glms),
                        Sample=1:(dim(inits)[3]))

init.means <- aaply(inits, c(1,2), mean)
init.sds <- aaply(inits, c(1,2), sd)
init.poptotals <- aaply(inits, c(2,3), sum)
### Generate Initial Conditions ###

# Using `spBayes`, I fit a model assuming constance average class distributions
# across the landscape but spatially autocorrelaation.  Future versions will
# Compare this to models with directional tendencies and possibly different
# forms of spatial correlation.



#plot.model <- InitModel(plot.sum)

# The function InitDraw uses the result of InitModel to draw a random initial
# populations at each point in the hex grid using a poisson distribution.  With
# error=TRUE, this also randomizes the expected value based on its standard
# erros.  The "array" argument will output an array suitable for initializing
# the population model, setting it FALSE gives us an easily plottable data
# frame
#init.test <- InitDraw(plot.model, error = FALSE, array=FALSE)

# Here I plot one 
# PLOTSPECIES <-"OTH"
# PLOTSIZECLASS <- "1"
# PLOTINFECT <- "S"
# 
# model.data <- subset(as.data.frame(init.test), 
#                      Species==PLOTSPECIES & SizeClass==PLOTSIZECLASS & 
#                        Infected==PLOTINFECT)
# real.data <- subset(plot.sum, 
#                     Species==PLOTSPECIES & SizeClass==PLOTSIZECLASS & 
#                       Infected==PLOTINFECT)
# poprange <- c(0, round_any(max(model.data$expected), 25, ceiling))
# radius = sqrt(AREA/pi)
# plot <- ggplot(real.data) + 
#   geom_hex(data=model.data, mapping=aes(x=Easting, y=Northing, fill=SimCount, 
#                                col=SimCount), 
#            stat="identity", size=0.1) + 
#   scale_x_continuous(limits=c(min(real.data$Easting) - radius, 
#                               max(real.data$Easting) + radius)) +
#   scale_y_continuous(limits=c(min(real.data$Northing) - radius, 
#                               max(real.data$Northing) + radius)) +
#   coord_equal() +
#   scale_color_gradient2(name = "Tree Count per Plot", limits=poprange) +
#   scale_fill_gradient2(guide="none", limits=poprange) +
#   theme(panel.grid=element_blank(), panel.background=element_blank())
# 
# for (i in 1:nrow(real.data)) {
#      plot <- plot + annotation_custom(grob=circleGrob(r=unit(0.5,"npc"), 
#                                                       gp=gpar(alpha=0.5, 
#                                                               fill="black")), 
#                                       xmin=real.data$Easting[i] - radius, 
#                                       xmax=real.data$Easting[i] + radius,
#                                       ymin=real.data$Northing[i] - radius,
#                                       ymax=real.data$Northing[i] + radius)
#      }
# plot <- plot + geom_text(mapping=aes(x=Easting, y=Northing, label=Count),
#                          data=real.data, col="white", cex=4)
# plot

## I need to improve the prediction on this using spBayes or geoRglm.  In the
# meantime, now we simulate the model using the output of InitDraw

#init <- InitDraw(plot.model, error = FALSE, array=TRUE)/AREA

## Let's initialize the model with the rest of the parameters

data(parms.Cobb2012.norm)
treeparms.df <- parms.Cobb2012.norm
treeparms.df <- within(treeparms.df, {
# Calculate the relative space requirements of tanoak size classes
# based on initial conditions
compete[1:4] <- 0.25*(sum(S.init[1:4])/S.init[1:4])

# Set recruitment rates to steady-state levels.  For Redwood and Bay, this is 
# simply the mortality rate divided by the density-dependence coefficient at 
# simulation start.  For tanoak, which has multiple size classes, it's a but
# more involved

S.recruit[5] <- S.mortality[5]/(1-sum(S.init * compete))
I.recruit[5] <- I.mortality[5]/(1-sum(S.init * compete))
S.recruit[6] <- S.mortality[6]/(1-sum(S.init * compete))
I.recruit[6] <- I.mortality[6]/(1-sum(S.init * compete))

A2 <- S.transition[1]/(S.transition[2] + S.mortality[2])
A3 <- (S.transition[2]/(S.transition[3] + S.mortality[3])) * A2
A4 <- (S.transition[3]/S.mortality[4]) * A3
S.recruit[1] <- (S.transition[1] + S.mortality[1])/(1-sum(S.init * compete)) -
                (S.recruit[2] * A2 + S.recruit[3] * A3 + S.recruit[4] * A4)


I.recruit[1] <- S.recruit[1] 
A2 <- NULL
A3 <- NULL
A4 <- NULL
})

## Get the locations matrix from plot.model

locations <- hex.pts

## Set time steps

time.steps <- 1:100

## Scale the dispersal kernel for the distance between plots

treeparms.df$kernel.fn <- "adjacent.dispersal"
treeparms.df$kernel.par1 <- 30
treeparms.df$kernel.par2 <- 33
treeparms.df$kernel.par3 <- 1
treeparms.df$kernel.par4 <- 0.5

init <- inits[,,1]

Rprof("SOD.out", line.profiling=TRUE, interval=0.02)
pop.df <- SODModel(treeparms.df,locations,time.steps,init, stochastic.d=FALSE,
                   df.out=FALSE, verbose=TRUE)
Rprof(NULL)

#inits <- foreach(i=1:4, .combine=abind2, .multicombine=TRUE) %do% {
#  InitDraw(plot.model, error = FALSE, array=TRUE)/AREA
#}

pops <- SODSims(inits=inits, parms.df=treeparms.df,locations=locations,time.steps=time.steps,stochastic.d=FALSE, .parallel=FALSE, verbose=TRUE)

out.to.df <- function(output) {
  out.df <- melt(output, value.name="Population")
  Classes <- as.character(out.df$Class)
  Classes <- matrix(unlist(strsplit(Classes, "[.,_]")), ncol=3,byrow=TRUE)
  colnames(Classes) <- c("Species", "SizeClass", "Infected")
  out.df$Class <- NULL
  out.df <- cbind(Classes, out.df)
  
}

pops.df <- out.to.df(pops)

pops.df$Key <- do.call(paste, c(pops.df[,1:6], sep="."))
  
require(data.table)

pops.dt <- data.table(pops.df, key="Key")

## Functions to summarize data

# Output array to data frame

# - Aggregate classes 

# - Large and small tanoak over time, plot each with error bars

# - Animate Given CLASS or 

#pop.df.totals <- ddply(pop.df, c("Time", "Disease", "Species", "SizeClass"),
#                       summarise, TotPop=mean(Population))
#dynamic.plot <- ggplot(pop.df.totals, 
#                       aes(x=Time,y=TotPop, fill=Disease, color=SizeClass)) + 
#                  geom_area(position="stack", alpha=0.6) + facet_grid(~Species)
#dynamic.plot