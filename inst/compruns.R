#Pseudocode for runs tomorrow.

# 1. Reproduce Cobb2012 infection results again
library(SODDr)
library(plyr)
library(stringr)
library(ggplot2)
library(grid)
library(hexbin)
data(parms.Cobb2012)
locations <- MakeLattice(nx=20,ny=20,dist=1)
initial.vec = as.vector(rbind(parms.Cobb2012$S.init, parms.Cobb2012$I.init))
init.Cobb2012 <- matrix(data=initial.vec, nrow=nrow(locations), 
               ncol=2*nrow(parms.Cobb2012), byrow=TRUE)
time.steps <- 1:100
init.Cobb2012[190,2] <- init.Cobb2012[190,1]
init.Cobb2012[190,1] <- 0
init.Cobb2012[190,10] <- init.Cobb2012[190,9]
init.Cobb2012[190,9] <- 0

pop <- SODModel(parms=parms.Cobb2012,locations=locations,times=time.steps,
                init=init.Cobb2012)
Cobb.2012.df <- merge(data.frame(Location=rownames(locations),locations),
                      melt(pop))
Cobb.2012.df$Size <- factor(ifelse(Cobb.2012.df$SizeClass %in% c("1","2"), 
                                   "Small", "Big"), 
                  levels=c("Small", "Big"))
df <- ddply(Cobb.2012.df, .(Replicate, Time, Species, Size, SizeClass, Infected),
            summarize, TotPop = sum(Population))
df <- ddply(df, .(Replicate, Time), transform, FracPop = TotPop/sum(TotPop))
df <- ddply(df, .(Replicate, Time, Species, Size), summarize, Frac=sum(FracPop))
plot.Cobb2012.LIDEsize <- ggplot(subset(df, Replicate==1 & Species=="LIDE")) +
  geom_line(mapping=aes(x=Time, y=Frac, col=Size), lwd=1.5) +
  labs(x="Time (years)", y="Fraction of Population") +
  scale_y_continuous(limit=c(0,0.6)) +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(colour="#ECECEC", size=0.5, linetype=1),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
        legend.text=element_text(size=22),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.Cobb2012.LIDEsize

times <- c(1,2,5,10,20,50)
df <- subset(Cobb.2012.df, Time %in% times & Species %in% c("LIDE", "UMCA"))
df <- ddply(df, .(Replicate, Time, Location, Infected), summarize, Tot = sum(Population))
df <- subset(ddply(df, .(Replicate, Time, Location), transform,
                   FracInf=Tot/sum(Tot)),
             Infected=="I")
df <- merge(data.frame(Location=rownames(locations),locations), df)

plot.Cobb2012.spread <- ggplot(df, aes(x=x, y=y, fill=FracInf)) +
  geom_tile() +
  facet_wrap(~Time) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_gradient(low="white", high="red", limits=c(0,1)) +
  coord_equal() +
  labs(x="x", y="y") + 
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
        #legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
       legend.text=element_text(size=18),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.Cobb2012.spread

# Do the same on the hexagonal grid and a larger carrying capacity
SUBSITE="SDC"
data(SOD.plots)
plot.sum <- plot.sums(SOD.plots, SUBSITE)
AREA=500
locations.hex <- hex.grid(plot.sum, hex.area=AREA)
init.hex <- matrix(50*init.Cobb2012[1,], ncol=ncol(init.Cobb2012),
                   nrow=nrow(locations.hex), byrow=TRUE)
init.hex[95,2] <- init.hex[95,1]
init.hex[95,1] <- 0
init.hex[95,10] <- init.hex[95,9]
init.hex[95,9] <- 0
#Change kernel
parms.hex <- parms.Cobb2012
parms.hex$kernel.par1 <- 5
parms.hex$kernel.par2 <- 30

pop.hex <- SODModel(parms=parms.hex,locations=locations.hex,times=time.steps,
                init=init.hex, K=50)

hex.df <- merge(data.frame(Location=rownames(locations.hex),locations.hex),
                      melt(pop.hex))
hex.df$Size <- factor(ifelse(hex.df$SizeClass %in% c("1","2"), 
                                   "Small", "Big"), 
                  levels=c("Small", "Big"))
df <- ddply(hex.df, .(Replicate, Time, Species, Size, SizeClass, Infected),
            summarize, MeanPop = mean(Population))
df <- ddply(df, .(Replicate, Time), transform, FracPop = MeanPop/sum(MeanPop))
df <- ddply(df, .(Replicate, Time, Species, Size), summarize, Frac=sum(FracPop))
plot.hex.LIDEsize <- ggplot(subset(df, Replicate==1 & Species=="LIDE")) +
  geom_line(mapping=aes(x=Time, y=Frac, col=Size), lwd=1.5) +
  labs(x="Time (years)", y="Fraction of Population") +
  scale_y_continuous(limit=c(0,0.6)) +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(colour="#ECECEC", size=0.5, linetype=1),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
        legend.text=element_text(size=22),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.hex.LIDEsize

times <- c(1,2,5,10,20,50)
df <- subset(hex.df, Time %in% times & Species %in% c("LIDE", "UMCA"))
df <- ddply(df, .(Replicate, Time, Location, Infected), summarize, Tot = sum(Population))
df <- subset(ddply(df, .(Replicate, Time, Location), transform,
                   FracInf=Tot/sum(Tot)),
             Infected=="I")
df <- merge(data.frame(Location=as.numeric(rownames(locations.hex)),
                                locations.hex),
            df)

plot.hex.spread <- ggplot(df, aes(x=Easting, y=Northing, col=FracInf, 
                                  fill=FracInf)) +
  geom_hex(stat="identity") +
  facet_wrap(~Time) +
  scale_fill_gradient(low="white", high="red", limits=c(0,1)) +
  scale_color_gradient(low="white", high="red", limits=c(0,1)) +
  coord_equal() +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
        #legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
       legend.text=element_text(size=18),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.hex.spread

## Do the same with the actual start values
require(spBayes)
SPM <- readRDS("SPMlong.rds")

## Record which classes have zero count and eliminate these
zero.classes <- unlist(dlply(plot.sum, .(Species, SizeClass, Infected), 
                                function(x) {sum(x$Count)==0}))
class.names <- names(zero.classes)
plot.sum2 <- subset(plot.sum, !(paste(Species, SizeClass, Infected, sep=".") == 
                                  names(zero.classes)[which(zero.classes)]))

## Create a hexagonal grid over the area of the data

locs <- hex.grid(plot.sum, AREA)

## Record essential variables
n.class <- length(class.names)
n.class2 <- sum(!zero.classes)
n.loc <- nrow(locations.hex)
n.plot <- length(unique(plot.sum$Plot))

#Set the number of final samples to draw:
n.inits <- 1000
n.batch <- 1000
batch.length <- 150
n.samples <- n.batch*batch.length
burn.in <- ceiling(0.9*n.samples)
sub.samps <- (burn.in+1):n.samples
spPredicts <- spPredict(SPM, pred.coords=locations.hex, 
                        pred.covars=mkMvX(rep(list(matrix(1, n.loc)), n.class2)),
                        start=burn.in, end=n.samples,
                        thin=floor((n.samples-burn.in)/n.inits),
                        verbose=TRUE)



#Reduce to n.inits by random sampling
sPy <- spPredicts$p.y.predictive.samples[, 
       sample(ncol(spPredicts$p.y.predictive.samples), n.inits)]

#Reorganize the data into
sPredicts.m <- aaply(1:n.loc, 1, function(x) {
  sPy[((x-1)*n.class2 + 1):(x*n.class2),]})
#sPredicts.m <- aperm(sPredicts.m, c(2,1,3))
n.z <- sum(zero.classes)
z.p <- array(0, dim=c(n.loc, n.z, n.inits))
sPredicts.m <- abind(sPredicts.m, z.p, along=2)
dimnames(sPredicts.m) <- list(Location=1:n.loc,
                              Class=c(class.names[which(!zero.classes)],
                                      class.names[which(zero.classes)]),
                              Replicate=1:n.inits)
sPredicts.m <- sPredicts.m[,class.names,]

#Generate inital conditions via 
sp.inits <- array(rpois(length(sPredicts.m), sPredicts.m), dim=dim(sPredicts.m),
                  dimnames=dimnames(sPredicts.m))
sp.stats <- abind(aaply(sPredicts.m, c(1,2), mean),
                  aaply(sPredicts.m, c(1,2), sd),
                  aaply(sp.inits, c(1,2), mean), 
                  aaply(sp.inits, c(1,2), sd),
                  aaply(sp.inits, c(1,2), max), 
                  aaply(sPredicts.m, c(1,2), median),
                  aaply(sp.inits, c(1,2), median),
                  aaply(sPredicts.m, c(1,2), function(x) {
                          fit <- density(x)
                          fit$x[which.max(fit$y)]
                        }),
                  along=3,
                  new.names=c(dimnames(sp.inits)[1:2], list(Stat=c("model.mean",
                                                         "model.sd",
                                                         "pred.mean",
                                                         "pred.sd",
                                                          "pred.max",
                                                          "model.med",
                                                          "pred.med",
                                                          "model.mode"))))
pop.sp <- SODModel(parms=parms.hex,locations=locations.hex,times=time.steps,
                init=sp.stats[,,"model.med"], K=50)

hex.sp <- merge(data.frame(Location=rownames(locations.hex),locations.hex),
                      melt(pop.sp))
hex.sp$Size <- factor(ifelse(hex.sp$SizeClass %in% c("1","2"), 
                                   "Small", "Big"), 
                  levels=c("Small", "Big"))
df <- ddply(hex.sp, .(Replicate, Time, Species, Size, SizeClass, Infected),
            summarize, MeanPop = mean(Population))
df <- ddply(df, .(Replicate, Time), transform, FracPop = MeanPop/sum(MeanPop))
df <- ddply(df, .(Replicate, Time, Species, Size), summarize, Frac=sum(FracPop))
plot.sp.LIDEsize <- ggplot(subset(df, Replicate==1 & Species=="LIDE")) +
  geom_line(mapping=aes(x=Time, y=Frac, col=Size), lwd=1.5) +
  labs(x="Time (years)", y="Fraction of Population") +
  scale_y_continuous(limit=c(0,0.6)) +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(colour="#ECECEC", size=0.5, linetype=1),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
        legend.text=element_text(size=22),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.sp.LIDEsize

times <- c(1,2,5,10,20,50)
df <- subset(hex.sp, Time %in% times & Species %in% c("LIDE", "UMCA"))
df <- ddply(df, .(Replicate, Time, Location, Infected), summarize, Tot = sum(Population))
df <- subset(ddply(df, .(Replicate, Time, Location), transform,
                   FracInf=Tot/sum(Tot)),
             Infected=="I")
df <- merge(data.frame(Location=as.numeric(rownames(locations.hex)),
                                locations.hex),
            df)

plot.sp.spread <- ggplot(df, aes(x=Easting, y=Northing, col=FracInf, 
                                  fill=FracInf)) +
  geom_hex(stat="identity") +
  facet_wrap(~Time) +
  scale_fill_gradient(low="white", high="red", limits=c(0,1)) +
  scale_color_gradient(low="white", high="red", limits=c(0,1)) +
  coord_equal() +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
        #legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
       legend.text=element_text(size=18),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.sp.spread

df <- subset(hex.sp, Time==1)
plot.sp.init <- ggplot(df, aes(x=Easting, y=Northing, col=Population,
                               fill=Population)) +
  geom_hex(stat="identity") +
  facet_wrap(~Species + SizeClass + Infected) +
  scale_color_gradient(trans="log") +
  scale_fill_gradient(trans="log") +
  coord_equal() +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
#        legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
        legend.text=element_text(size=18),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.sp.init

df <- subset(hex.sp, Time==1)
df <- subset(ddply(df, .(Location, Species, SizeClass), transform, 
                   PctInf=Population/sum(Population)), Infected=="I")
plot.sp.initInf <- ggplot(df, aes(x=Easting, y=Northing, col=PctInf,
                               fill=PctInf)) +
  geom_hex(stat="identity") +
  facet_wrap(~Species + SizeClass) +
  scale_color_gradient(low="white", high="red", limits=c(0,1)) +
  scale_fill_gradient(low="white", high="red", limits=c(0,1)) +
  coord_equal() +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.title=element_blank(),
#        legend.position=c(0.75,0.6),
        legend.key=element_rect(fill="white"),
        legend.key.size=unit(1.5, "cm"),
        legend.text=element_text(size=18),
        axis.title=element_text(size=24),
#        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))
plot.sp.initInf

packs <- .packages()
save.image(file="cmpruns.Rdata")

### Transfer to EC2

load("cmpruns.Rdata")
lapply(packs,function(x){library(x,character.only=TRUE)}) 

pop.hexsch <- SODModel(parms=parms.hex,locations=locations.hex,times=time.steps,
                       init=init.hex, stochastic.d=TRUE, reps=50, K=50)