library(SODDr)
library(hexbin)
library(spBayes)
library(plyr)
library(reshape2)
library(grid)
library(ggplot2)

# Set some essential variables
SUBSITE <- "SDC"  # Which subsite of the plot data 
AREA <- 500

# Load the data
data(SOD.plots)

# Summarize data by subsite
plot.sum <- plot.sums(SOD.plots, subsite=SUBSITE, year="2002")

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
n.loc <- nrow(locs)
n.plot <- length(unique(plot.sum$Plot))

## Fit GLMs to the data in each class to get initial estimates
plot.glms <- dlply(plot.sum2, .(Species, SizeClass, Infected), function(x) {
                     glm(Count ~ 1, family="poisson", data=x)
                    })

n.class2 <- length(plot.glms)


## Set MCMC parameters
n.batch <- 1000
batch.length <- 150
n.samples <- n.batch*batch.length
burn.in <- ceiling(0.9*n.samples)
sub.samps <- (burn.in+1):n.samples


#Set up modeling structure
formulas <- alply(paste0(names(plot.glms),"~1"),1,as.formula)
betas <- unlist(lapply(plot.glms, coef))
names(betas) <- NULL
phis <- rep(0.01, n.class2)
As <- chol(diag(unlist(lapply(plot.glms,
                              vcov))))[lower.tri(diag(1,n.class2), TRUE)]
            
##Restrucutre data for use with spBayes
plot.cast <- dcast(plot.sum, Plot + Site + Subsite + Easting + Northing ~ 
                     Species + SizeClass + Infected, value.var="Count")
names(plot.cast) <- gsub("_",".", names(plot.cast))

##Load previous file:

start.file <- "SPM.rds"
if (file.exists(start.file)) {
  start.data <- restart.mcmc(readRDS(start.file))
  starting=start.data$starting
  tuning=start.data$tuning
} else {
  starting=list("beta"=betas,"phi"=phis,"A"=As, "w"=0)
  tuning=list("beta"=rep(0.5,n.class2), "phi"=rep(4,n.class2),
                           "A"=rep(4, length(As)), "w"=0.1)
}

## Fit the model via MCMC
SPMlong <- spMvGLM(formulas, family="poisson", data=plot.cast,
               coords=as.matrix(plot.cast[,c("Easting","Northing")]),
               starting=starting, tuning=tuning,
               priors=list("beta.Flat", "phi.Unif"=list(rep(0.00001,n.class2), 
                                                        rep(0.6,n.class2)),
                           "K.IW"=list(n.class2+1, diag(0.1,n.class2))),
               cov.model="exponential",
               amcmc=list("n.batch"=n.batch, "batch.length"=batch.length,
                          "accept.rate"=0.43),
               verbose=TRUE, n.report=1)

#A function to draw MCMC initial conditions from the previous run
saveRDS(SPMlong, "SPMlong.rds")

SPM <- readRDS("SPMlong.rds")

#Set the number of final samples to draw:
n.inits <- 1000

#Generate predictions
spPredicts <- spPredict(SPM, pred.coords=locs, 
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
                  aaply(sp.inits, c(1,2), median), along=3,
                  new.names=c(dimnames(sp.inits)[1:2], list(Stat=c("model.mean",
                                                         "model.sd",
                                                         "pred.mean",
                                                         "pred.sd",
                                                          "pred.max",
                                                          "model.med",
                                                          "pred.med"))))
sp.stats <- melt(sp.stats, value.name="Population")
colnames(sp.stats)[1:3] <- c(names(dimnames(sp.inits))[1:2], "Stat")
sp.stats <- merge(data.frame(Location=rownames(locs), locs), sp.stats)

VIZ.SP <- "OTH"
VIZ.SIZE <- 1
VIZ.INF <- "S"
plotclass <- paste(VIZ.SP, VIZ.SIZE, VIZ.INF, sep=".")
plots.df <- subset(plot.sum, 
                   Species==VIZ.SP & SizeClass==VIZ.SIZE & Infected==VIZ.INF)
VIZ.STAT <-"pred.med"
radius <- sqrt(AREA/pi)
meanplot <- ggplot(data=sp.stats, mapping=aes(x=Easting, y=Northing)) +
  geom_hex(data=subset(sp.stats, Class==plotclass & Stat==VIZ.STAT),
           stat="identity",
           mapping=aes(x=Easting, y=Northing, fill=Population, col=Population)) +
  scale_x_continuous(limits=c(min(plots.df$Easting) - radius, 
                              max(plots.df$Easting) + radius),
                     name="UTM Easting (m)") +
  scale_y_continuous(limits=c(min(plots.df$Northing) - radius, 
                              max(plots.df$Northing) + radius),
                     name="UTM Northing (m)") +
  scale_fill_gradientn(colours=c("white","blue"), space="rgb", name=VIZ.STAT) +
  scale_color_gradientn(colours=c("white","blue"), space="rgb", name=VIZ.STAT) +
  coord_equal() +
  theme(panel.grid=element_blank(), panel.background=element_blank())

for (i in 1:nrow(plots.df)) {
     meanplot <- meanplot + annotation_custom(grob=circleGrob(r=unit(0.5,"npc"), 
                                                      gp=gpar(alpha=0.8,
                                                              fill="black")),
                                      xmin=plots.df$Easting[i] - radius, 
                                      xmax=plots.df$Easting[i] + radius,
                                      ymin=plots.df$Northing[i] - radius,
                                      ymax=plots.df$Northing[i] + radius)
     }

meanplot <- meanplot + geom_text(data=plots.df, 
                         mapping=aes(x=Easting, y=Northing, label=Count), 
                         col="white", cex=4)
meanplot



##### Parameter density plots

matplot(t(SPM$acceptance), type="l")


SPMout <- SPM$p.beta.theta.samples
Alocs <- grep("K\\[", dimnames(SPMout)[[2]])
philocs <- grep("phi\\[", dimnames(SPMout)[[2]])
betalocs <- grep("\\(Intercept\\)", dimnames(SPMout)[[2]])
dimnames(SPMout)[[2]][1:n.class2] <- paste0("beta.",names(plot.glms))
class.means <- melt(SPMout[sub.samps,betalocs])
names(class.means) <- c("Sample", "Class", "Parameter")
parmplot <- ggplot(class.means, aes(x=exp(Parameter), fill=Class)) +
  geom_density(alpha=0.8, col=NA) +
  facet_wrap(~Class, scales="free_y") +
  expand_limits(y=0) + 
  scale_x_continuous(name="Log of Mean Number of Trees per Cell") + 
  scale_y_continuous(name="Posterior Density") + 
  theme(panel.grid=element_blank(), panel.background=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank(),
        panel.border=element_rect(fill=NA, colour="grey"),
        legend.position="none")
parmplot

plot(SPM$p.beta.theta.samples)


acc.out <- as.data.frame(SPM$acceptance)
colnames(acc.out) <- 1:(ncol(SPM$acceptance))
acc.out$VAR <- dimnames(SPMout)[[2]]
acc.out$Group <- NA
acc.out$Group[betalocs] <- "betas"
acc.out$Group[Alocs] <- "As"
acc.out$Group[philocs] <- "phi"
acc.out$Group <- as.factor(acc.out$Group)
acc.out <- melt(acc.out)
acc.out$variable = as.numeric(acc.out$variable)
names(acc.out)[3:4] = c("Batch", "Acceptance")

acceptplot <- ggplot(acc.out, aes(x=Batch, y=Acceptance, col=Group, group=VAR)) +
  geom_line() + facet_wrap(~Group, nrow=3)
acceptplot