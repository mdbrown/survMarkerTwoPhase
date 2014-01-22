library(survMarkerTwoPhase)

source("simulateData.R")
Nsim = 1000
## semi parametric
ncc.sp.sim <- replicate( n = Nsim, 
                         run_one_sim(sample.design = "ncc", 
                                     estimation.method = "SP"))

save(ncc.sp.sim, file = "ncc.sp.sim.Rdata")
