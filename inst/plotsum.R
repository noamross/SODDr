require(plyr)
require(ggplot2)
require(sp)
require(hexbin)
require(grid)
require(SODDr)
require(compiler)

enableJIT(3)

# Set some essential variables
SUBSITE <- "Old Growth"  # Which subsite of the plot data 
AREA <- 500

# Load the data and merge subsites names
rawdata <- read.csv("data/redwood_plot_network_full.csv", comment.char="#")
subsites <- read.csv("data/plot-site list.csv")
names(subsites)[2] <- "Subsite"
rawdata <- merge(rawdata, subsites)

### Data Conversion ###

# Get initial census of live trees at the chosen subsite
census.2002 <- subset(rawdata, Year==2002 & Mortality==0 & Subsite==SUBSITE)

# Aggregate non-sporulating species and assume no infection in these
sporulators <- c("UMCA", "LIDE")
census.2002$Species <- factor(census.2002$Species, 
                              levels = c(levels(census.2002$Species), "OTH"))
census.2002$Species[which(!(census.2002$Species %in% sporulators))] <- "OTH"
census.2002$Infected[which(!(census.2002$Species %in% sporulators))] <- FALSE
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

### Generate Initial Conditions ###

# InitModel tries out a variety of spaital models to fit the tree distributions
# of each class and returns a SpatialPointsDataFrame describing a hex grid
# of expected populations in the plot.  The models are stored as attributes
# of the output.  Note that it throws a bunch of errors when it tries models
# that fail
plot.model <- InitModel(plot.sum)

# The function InitDraw uses the result of InitModel to draw a random initial
# populations at each point in the hex grid using a poisson distribution.  With
# error=TRUE, this also randomizes the expected value based on its standard
# erros.  The "array" argument will output an array suitable for initializing
# the population model, setting it FALSE gives us an easily plottable data
# frame
init.test <- InitDraw(plot.model, error = FALSE, array=FALSE)

# Here I plot one 
PLOTSPECIES <-"LIDE"
PLOTSIZECLASS <- "3"
PLOTINFECT <- FALSE

model.data <- subset(as.data.frame(init.test), 
                     Species==PLOTSPECIES & SizeClass==PLOTSIZECLASS & 
                       Infected==PLOTINFECT)
real.data <- subset(plot.sum, 
                    Species==PLOTSPECIES & SizeClass==PLOTSIZECLASS & 
                      Infected==PLOTINFECT)
poprange <- c(0, round_any(max(model.data$expected), 25, ceiling))
plot <- ggplot(real.data) + 
  geom_hex(data=model.data, mapping=aes(x=Easting, y=Northing, fill=SimCount, 
                               col=SimCount), 
           stat="identity", size=0.1) + 
  scale_x_continuous(limits=c(min(real.data$Easting) - radius, 
                              max(real.data$Easting) + radius)) +
  scale_y_continuous(limits=c(min(real.data$Northing) - radius, 
                              max(real.data$Northing) + radius)) +
  coord_equal() +
  scale_color_gradient2(name = "Tree Count per Plot", limits=poprange) +
  scale_fill_gradient2(guide="none", limits=poprange) +
  theme(panel.grid=element_blank(), panel.background=element_blank())

for (i in 1:nrow(real.data)) {
     plot <- plot + annotation_custom(grob=circleGrob(r=unit(0.5,"npc"), 
                                                      gp=gpar(alpha=0.5, 
                                                              fill="black")), 
                                      xmin=real.data$Easting[i] - radius, 
                                      xmax=real.data$Easting[i] + radius,
                                      ymin=real.data$Northing[i] - radius,
                                      ymax=real.data$Northing[i] + radius)
     }
plot <- plot + geom_text(mapping=aes(x=Easting, y=Northing, label=Count),
                         data=real.data, col="white", cex=4)
plot

## I need to improve the prediction on this using spBayes or geoRglm.  In the
# meantime, now we simulate the model using the output of InitDraw

init <- InitDraw(plot.model, error = FALSE, array=TRUE)/AREA

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

locations <- attr(plot.model, "locations")

## Set time steps

time.steps <- 1:100

## Scale the dispersal kernel for the distance between plots

treeparms.df$kernel.fn <- "adjacent.dispersal"
treeparms.df$kernel.par1 <- 30
treeparms.df$kernel.par2 <- 33
treeparms.df$kernel.par3 <- 1
treeparms.df$kernel.par4 <- 0.5

Rprof("SOD.out", line.profiling=TRUE, interval=0.02)
pop.df <- SODModel(treeparms.df,locations,time.steps,init, stochastic.d=FALSE)
Rprof(NULL)

pop.df.totals <- ddply(pop.df, c("Time", "Disease", "Species", "SizeClass"),
                       summarise, TotPop=mean(Population))
dynamic.plot <- ggplot(pop.df.totals, 
                       aes(x=Time,y=TotPop, fill=Disease, color=SizeClass)) + 
                  geom_area(position="stack", alpha=0.6) + facet_grid(~Species)
dynamic.plot