require(compiler)
require(gdata)
require(plyr)
require(reshape2)
require(ggplot2)

# Load all the function files
for (file in list.files("../R/", full.names=TRUE)) source(file)



#1.  Input and calculate the species parameters
treeparms.df <- read.csv("paper_tree_parms_eq.csv", stringsAsFactors=FALSE)

# Calculate missing values

# Calculate the relative space requirements of tanoak
treeparms.df$space[1:4] <- 0.25*.199*1.5/(treeparms.df$S.init[1:4])

# Set recruitment rates to steady-state levels
treeparms.df$S.recruit[5] <- treeparms.df$S.mortality[5]/(1-sum(treeparms.df$S.init * treeparms.df$space))
treeparms.df$I.recruit[5] <- treeparms.df$I.mortality[5]/(1-sum(treeparms.df$S.init * treeparms.df$space))
treeparms.df$S.recruit[6] <- treeparms.df$S.mortality[6]/(1-sum(treeparms.df$S.init * treeparms.df$space))
treeparms.df$I.recruit[6] <- treeparms.df$I.mortality[6]/(1-sum(treeparms.df$S.init * treeparms.df$space))

A2 <- treeparms.df$S.transition[1]/(treeparms.df$S.transition[2] + treeparms.df$S.mortality[2])
A3 <- A2 * treeparms.df$S.transition[2]/(treeparms.df$S.transition[3] + treeparms.df$S.mortality[3])
A4 <- A3* treeparms.df$S.transition[3]/treeparms.df$S.mortality[4]
treeparms.df$S.recruit[1] <- (treeparms.df$S.transition[1] + treeparms.df$S.mortality[1])/(1-sum(treeparms.df$S.init * treeparms.df$space)) -
  (treeparms.df$S.recruit[2] * A2 + treeparms.df$S.recruit[3] * A3 + treeparms.df$S.recruit[4] * A4)

treeparms.df$I.recruit[1] <- treeparms.df$S.recruit[1] 

# This dispersal kernel defined in this matrix is as follows:

#'This function returns the first argument for distances < 0.5, the second for distances [0.5,1.5], zero for others
adjacent.dispersal <- function(distance, local, adjacent) {
  ifelse(distance < 0.5, local,
         ifelse(distance < 1.5, adjacent,0)
  )
}


locations <- MakeLattice(20,20,1)

initial.vec = as.vector(rbind(treeparms.df$S.init, treeparms.df$I.init))
init <- matrix(data=initial.vec, nrow=nrow(locations), ncol=2*nrow(treeparms.df), byrow=TRUE)

#Infect trees in one location
#pop.init[190,2] <- pop.init[190,1]
#pop.init[190,1] <- 0
#pop.init[190,10] <- pop.init[190,9]
#pop.init[190,9] <- 0

time.steps <- 1:100

pop.df <- SODModel(treeparms.df,locations,time.steps,init)
#Rprof(NULL)

# Plot

# Convert to a data frame


densfun <- function(x) {sum(x)/nrow(locations)}
pop.df.totals <- ddply(pop.df, c("Time", "Disease","Species","AgeClass"), summarise, AvePop=sum(Population)/nrow(locations))
dynamic.plot <- ggplot(pop.df.totals, aes(x=Time,y=AvePop, fill=Disease, color=AgeClass)) + geom_area(position="stack", alpha=0.6) + facet_grid(~Species)
dynamic.plot
final.densities <-ddply(subset(pop.df, Time==tail(time.steps,1)), c("Species","AgeClass","Disease"), summarise, FinalPop=sum(Population)/nrow(locations))
paper.df <- ddply(transform(subset(pop.df.totals, Disease=="S"), Size = ifelse(AgeClass == "1" | AgeClass == "2", "Small", "Large")), c("Time", "Species", "Size"), summarise, SizePop=sum(AvePop))
paper.df <- ddply(paper.df, "Time", transform, Pct=SizePop/sum(SizePop))
paper.plot <- ggplot(subset(paper.df, Species==1), aes(x=Time, y=Pct, lty=Size)) + geom_line(lwd=1) + scale_y_continuous(limits=c(0,0.6), expand=c(0,0)) + scale_linetype_manual(values=c(1,5))
paper.plot
