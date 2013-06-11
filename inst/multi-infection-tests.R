#' First test of a model with multiple infection
#'TODO: Convert to continuous
#'TODO: Add density dependence.  
#'
require(reshape2)
require(ggplot2)
require(deSolve)
require(noamtools)

parms <- c( 
  f_j=0.01,
  f_a=0.01,
  g=0.1,
  d_j=0.005,
  d_a=0.005,
  alpha=0.05,
  lambda=0.3,
  K=50,
  mu=0.01
 # H_0 = 50
  )
  
basic.model <- function(t, y, parms) {
  list2env(as.list(y), environment())
  list2env(as.list(parms), environment())
  dJ <- A * f_a + J * (f_j - d_j - g) - alpha * PJ
  dA <- J * g - (A * d_a) - alpha * PA
  dPJ <- lambda * (PJ + PA) * J/(J+A) - PJ*(d_j + g + alpha*(1 + PJ/J))
  dPA <- lambda * (PJ + PA) * A/(J+A) + PJ * g - PA*(d_a + alpha*(1 + PA/A))
  return(list(c(dJ=dJ, dA=dA, dPJ=dPJ, dPA=dPA)))
}

dd.model <- function(t, y, parms) {
  list2env(as.list(y), environment())
  list2env(as.list(parms), environment())
  dJ <- A * f_a * (1 - (J+A)/K) + J * (f_j * (1 - (J+A)/K) - d_j - g) - alpha * PJ
  dA <- J * g - (A * d_a) - alpha * PA
  dPJ <- lambda * (PJ + PA) * J/(J+A) - PJ*(d_j + g + alpha*(1 + PJ/J))
  dPA <- lambda * (PJ + PA) * A/(J+A) + PJ * g - PA*(d_a + alpha*(1 + PA/A))
  return(list(c(dJ=dJ, dA=dA, dPJ=dPJ, dPA=dPA)))
}

ddmu.model <- function(t, y, parms) {
  list2env(as.list(y), environment())
  list2env(as.list(parms), environment())
  dJ <- A * f_a * (1 - (J+A)/K) + J * (f_j * (1 - (J+A)/K) - d_j - g) - alpha * PJ
  dA <- J * g - (A * d_a) - alpha * PA
  dPJ <- lambda * (PJ + PA) * J/K - PJ*(d_j + mu + g + alpha*(1 + PJ/J))
  dPA <- lambda * (PJ + PA) * A/K + PJ * g - PA*(d_a + mu + alpha*(1 + PA/A))
  return(list(c(dJ=dJ, dA=dA, dPJ=dPJ, dPA=dPA), c(dJ=dJ, dA=dA, dPJ=dPJ, dPA=dPA)))
}

init <- c(J=1.190476, A=23.8108, PJ=0.119, PA=2.381)

out <- as.data.frame(lsoda(y=init, times=1:100, func=ddmu.model, parms=parms))
df <- as.data.frame(out)
names(df)[1] <- "Time"

df <- within(df, {
  list2env(as.list(parms), environment())
  pctJ <- J/(J + A)
  pctA <- A/(J + A)
  nJ <- PJ / J
  nA <- PA / A
  J.inf <- 1 - exp(-nJ)
  A.inf <- 1 - exp(-nA)
  Inf.dens <- (J*J.inf + A*A.inf)*20
  J.mort <- d_j + alpha * PJ / (J * (1 - J.inf))
  A.mort <- d_a + alpha * PA / (A * (1 - A.inf))
  J.yrs <- 1/J.mort
  A.yrs <- 1/A.mort
  J.infrate <- 1 - exp(-lambda * (PJ + PA) * J * exp(-PJ/J) / K)
  A.infrate <- 1 - exp(-lambda * (PJ + PA) * A * exp(-PA/A) / K)
  J.infyrs <- 1/J.infrate
  A.infyrs <- 1/A.infrate
  rm(list=names(parms))
})

df <- df[,sort(names(df), decreasing=TRUE)]

mdf <- melt(df, id.vars=c("Time", "Inf.dens"), variable.name="Class", value.name="Population")
JAlab <- scale_color_discrete(labels=c("Small Trees","Big Trees"))
ggplot(subset(mdf, Class %in% c("J", "A")), aes(x=Time, y=Population, col=Class)) + geom_line(lwd=1) + theme_nr + ylab("Population") + JAlab
ggplot(subset(mdf, Class %in% c("pctJ", "pctA")), aes(x=Time, y=Population, col=Class)) + geom_line(lwd=1) + theme_nr + ylab("Fraction of Population") + JAlab
ggplot(subset(mdf, Class %in% c("PJ", "PA")), aes(x=Time, y=Population, col=Class)) + geom_line(lwd=1) + theme_nr + ylab("Number of Infections") + JAlab
ggplot(subset(mdf, Class %in% c("nJ", "nA")), aes(x=Time, y=Population, col=Class)) + geom_line(lwd=1) + theme_nr + ylab("Infections per Individual") + JAlab
ggplot(subset(mdf, Class %in% c("J.inf", "A.inf")), aes(x=Time, y=Population, col=Class)) + geom_line(lwd=1) + theme_nr + ylab("Fraction infected") + JAlab
ggplot(subset(mdf, Class %in% c("J.mort", "A.mort")), aes(x=Time, y=Population, col=Class)) + geom_line(lwd=1) + theme_nr + ylab("Mortality Rate of Infected Individuals") + JAlab

ggplot(subset(mdf, Class %in% c("J.yrs", "A.yrs")), aes(x=Inf.dens, y=Population, col=Class)) + geom_point(cex=4) + theme_nr + xlab("Density of infected trees (1/ha)") + ylab("Years to death of infected individuals") + JAlab
ggplot(subset(mdf, Class %in% c("J.infyrs", "A.infyrs") ), aes(x=Inf.dens, y=Population, col=Class)) + geom_point(cex=4) + theme_nr + xlab("Density of infected trees (1/ha)") + ylab("Years to infection") + JAlab