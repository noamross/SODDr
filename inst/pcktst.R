# package tests May 9, 2013
require(doMC)
registerDoMC()
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
init <- init*50
init2 <- init
init2[190,2] <- init2[190,1]
init2[190,1] <- 0
init2[190,10] <- init2[190,9]
init2[190,9] <- 0
inits <- abind2(init,init2)


pop <- SODModel(parms=parms.Cobb2012,locations=locs,times=time.steps,init=sp.med,
                   stochastic.d=TRUE, parallel=FALSE, K=50)
aaply(pop[,,,1], c(1,3), mean)
library(ggplot2)
require(grid)
library(plyr)
pop.df <- melt(pop)
pop.df$Size <- factor(ifelse(pop.df$SizeClass %in% c(1,2), "Small", "Big"), 
                  levels=c("Small", "Big"))
df <- ddply(pop.df, .(Replicate, Time, Species, Size, Infected), summarize, 
          MeanPop = mean(Population))
plot <- ggplot(subset(df, Replicate==1)) +
  geom_line(subset(df, Replicate==1), mapping=aes(x=Time, y=MeanPop, col=Species, lwd=Size, lty=Infected)) +
  labs(x="Time (years)", y="Mean Stem Count/Plot") +
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
 


##########

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
                                    SPG$p.beta.theta.samples[sub.samps,2])
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

spinits <- laply(sPredicts, function(SPR) {
                apply(SPR$p.y.predictive.sample,1, function(x) {
                        x <- sample(x, n.inits)
                        rpois(n.inits, x) 
                        })
                })
spinits <- aperm(spinits, c(3,1,2))
dimnames(spinits) <- list(Location=1:(dim(spinits)[1]),Class=names(spGLMs),
                        Sample=1:(dim(spinits)[3]))
summary(spinits[,,1])


sppops <- SODModel(init=spinits[,,1], parms=parms.Cobb2012,locations=locations,
                times=time.steps, stochastic.d=FALSE, parallel=FALSE)

df2 <- melt(sppops)
df2$Replicate <- 1
df2$Size <- factor(ifelse(df2$SizeClass %in% c(1,2), "Small", "Big"), 
                  levels=c("Small", "Big"))
df2 <- ddply(df2, .(Replicate, Time, Species, Size, Infected), summarize, 
          MeanPop = mean(Population))
plot <- ggplot(subset(df2, Replicate==1)) +
  geom_line(subset(df2, Replicate==1), mapping=aes(x=Time, y=MeanPop, col=Species, lwd=Size, lty=Infected)) +
  labs(x="Time (years)", y="Mean Stem Count/Plot") +
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
 
