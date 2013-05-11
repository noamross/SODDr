# package tests May 9, 2013

require(SODDr)
data(parms.Cobb2012)
data(SOD.plots)
plot.sum <- plot.sums(SOD.plots, "SDC")
AREA=500
locations <- hex.grid(plot.sum, hex.area=AREA)
rm(SOD.plots)
time.steps <- 1:100


# Set dispersal kernel parameters to match grid size
parms.Cobb2012$kernel.par1 <- 20
parms.Cobb2012$kernel.par1 <- 30

initial.vec = as.vector(rbind(parms.Cobb2012$S.init, parms.Cobb2012$I.init))
init <- matrix(data=initial.vec, nrow=nrow(locations), 
               ncol=2*nrow(parms.Cobb2012), byrow=TRUE)

pop <- SODModel(parms=parms.Cobb2012,locations=locations,times=time.steps,init=init,
                   stochastic.d=FALSE)
library(ggplot2)
library(plyr)
pop.df <- melt(pop)
pop.df.totals <- ddply(pop.df, c("Time", "Infected", "Species", "SizeClass"),
                       summarise, TotPop=mean(Population))
dynamic.plot <- ggplot(pop.df.totals, 
                       aes(x=Time,y=TotPop, fill=Infected, color=SizeClass)) + 
                  geom_area(position="stack", alpha=0.6) + facet_grid(~Species)
dynamic.plot
init2 <- init

init2[190,2] <- init2[190,1]
init2[190,1] <- 0
init2[190,10] <- init2[190,9]
init2[190,9] <- 0
inits <- abind2(init,init2)

pops2 <- SODModel(init=init2, reps=10, parms=parms.Cobb2012,locations=locations,
                times=time.steps, verbose=TRUE)
pops2.df <- melt(pops2)

pops2.df$Size <- factor(ifelse(pops2.df$SizeClass %in% c(1,2), "Small", "Big"), 
                  levels=c("Small", "Big"))
df <- ddply(pops2.df, .(Replicate, Time, Species, Size), summarize, 
            TotPop = sum(Population))
df <- ddply(df, .(Replicate, Time), summarize, 
            Species=Species, Size=Size, PctPop=TotPop/sum(TotPop))
df.means <- ddply(df, .(Time, Species, Size), summarize, MeanPop=mean(PctPop))
df$Run=factor(paste0(df$Replicate,".",df$Size))

plot <- ggplot(subset(df, Species=="LIDE")) +
  geom_line(data=subset(df, Species=="LIDE"), 
            mapping=aes(x=Time, y=PctPop, col=Size, group=Run), alpha=0.2) +
  geom_line(data=subset(df.means, Species=="LIDE"), 
            mapping=aes(x=Time, y=MeanPop, col=Size),  lwd=1) +
  labs(x="Time (years)", y="Tanoak Population") +
  scale_color_discrete(labels=c("Small Tanoak", "Large Tanoak")) +
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
plot

require(spBayes)
plot.glms <- dlply(plot.sum, .(Species, SizeClass, Infected), function(x) {
                     glm(Count ~ 1, family="poisson", data=x)
                    })

nclass <- length(plot.glms)
nplots <- length(unique(plot.sum$Plot))
n.batch <- 1500
batch.length <- 100
n.samples <- n.batch*batch.length
n.inits <- 10
burn.in <- ceiling(0.9*n.samples)
sub.samps <- (burn.in+1):n.samples

spGLMs <- llply(plot.glms, function(GLM) {
                  spGLM(Count~1, family="poisson", data=GLM$data,
                  coords=as.matrix(GLM$data[,c("Easting", "Northing")]),
        starting=list("beta"=coef(GLM),"phi"=0.001,"sigma.sq"=vcov(GLM), "w"=0),
        priors=list("beta.Flat", "phi.Unif"=c(0.000001, 0.03), "sigma.sq.IG"=c(2, 1)),
        tuning=list("beta"=0.1, "phi"=2, "sigma.sq"=0.5, "w"=0.5),
        amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
        cov.model="exponential", verbose=TRUE)
                })

class.means <- ldply(spGLMs, function(SPG) {
                       data.frame(MeanLogCount=
                                    SPG$p.beta.theta.samples[sub.samps,1])
                       })
names(class.means)[1] <- "Class"
parmplot <- ggplot(subset(class.means, Class!="OTH.1.I"), 
              aes(x=MeanLogCount, fill=Class)) +
  geom_density(alpha=0.8, col=NA) +
  facet_wrap(~Class, scales="free") +
  expand_limits(y=0) + 
  scale_x_continuous(name="Log of Mean Number of Trees per Cell") + 
  scale_y_continuous(name="Posterior Density") + 
  theme(panel.grid=element_blank(), panel.background=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank(),
        panel.border=element_rect(fill=NA, colour="grey"),
        legend.position="none")
parmplot

sPredicts <- llply(spGLMs, function(SPG) {
                     spPredict(SPG, pred.coords=locations, 
                               pred.covars=matrix(1, nrow(locations)),
                               start=burn.in, end=n.samples,
                               thin=floor((n.samples-burn.in)/n.inits),
                               verbose=FALSE)
                  })

inits <- laply(sPredicts, function(SPR) {
                apply(SPR$p.y.predictive.sample,1, function(x) {
                        x <- sample(x, n.inits)
 #                       rpois(n.inits, x) 
                        })
                })
inits <- aperm(inits, c(3,1,2))
dimnames(inits) <- list(Location=1:(dim(inits)[1]),Class=names(spGLMs),
                        Sample=1:(dim(inits)[3]))
summary(inits[,,4])

inits <- inits/AREA

pops3 <- SODModel(init=inits, parms=parms.Cobb2012,locations=locations,
                times=time.steps)
pops3.df <- melt(pops3)

pops3.df$Size <- factor(ifelse(pops3.df$SizeClass %in% c(1,2), "Small", "Big"), 
                  levels=c("Small", "Big"))
df <- ddply(pops3.df, .(Replicate, Time, Species, Size), summarize, 
            TotPop = sum(Population))
df <- ddply(df, .(Replicate, Time), summarize, 
            Species=Species, Size=Size, PctPop=TotPop/sum(TotPop))
df.means <- ddply(df, .(Time, Species, Size), summarize, MeanPop=mean(PctPop))
df$Run=factor(paste0(df$Replicate,".",df$Size))

plot <- ggplot(subset(df, Species=="LIDE")) +
  geom_line(data=subset(df, Species=="LIDE"), 
            mapping=aes(x=Time, y=PctPop, col=Size, group=Run), alpha=0.2) +
  geom_line(data=subset(df.means, Species=="LIDE"), 
            mapping=aes(x=Time, y=MeanPop, col=Size),  lwd=1) +
  labs(x="Time (years)", y="Tanoak Population") +
  scale_color_discrete(labels=c("Small Tanoak", "Large Tanoak")) +
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
plot
