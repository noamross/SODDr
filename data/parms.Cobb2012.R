# This script takes the data from parms.Cobb_2012_raw and calculates
# the competitive coefficients

parms.Cobb2012 <- utils::read.csv("parms.Cobb2012_raw.csv",
                                  stringsAsFactors=FALSE)
parms.Cobb2012 <- within(parms.Cobb2012, {
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