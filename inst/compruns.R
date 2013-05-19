#Pseudocode for runs tomorrow.

# 1. Reproduce Cobb2012 infection results again
library(SODDr)
data(parms.Cobb2012)
locations <- MakeLattice(nx=20,ny=20,dist=1)
initial.vec = as.vector(rbind(parms.Cobb2012$S.init, parms.Cobb2012$I.init))
init <- matrix(data=initial.vec, nrow=nrow(locations), 
               ncol=2*nrow(parms.Cobb2012), byrow=TRUE)
time.steps <- 1:100
init2 <- init
init2[190,2] <- init2[190,1]
init2[190,1] <- 0
init2[190,10] <- init2[190,9]
init2[190,9] <- 0
inits <- abind2(init,init2)

pop <- SODModel(parms=parms.Cobb2012,locations=locations,times=time.steps,init=inits,
                stochastic.d=FALSE, parallel=FALSE, K=50)

# 2. Reproduce Cobb2012 with 100 runs of demographic stochasticity

pop.sch <- SODModel(parms=parms.Cobb2012,locations=locations,times=time.steps,init=inits,
                       stochastic.d=TRUE, parallel=TRUE, K=50, reps=25)

saveRDS(pop.sch, "pop.sch.rds")


# 3. Reproduce Cobb2012 with the SDC hex grid

# 4. Reproduce Cobb2021 with the SDC hex grid and demographic stochasticity

# 5. Reproduce Cobb2012 with SDC hex grid using median spBayes predictions 
#  Main question:  Why instantaneous infection?
#   Possibility: infected UMCA already across the landscape?
# 6.  Reproduce Cobb 2012 with SDC hex grid and median spBayes predictions, demog stoch
# 7.  Reproduce Cobb 2012 with SDC hex grid and randomized spBayes predictions.  


## Coding: pull out initial work from loop
## Coding: make robust to run failure/error.  å