library(survMarkerTwoPhase)

source("simulateData.R")
Nsim = 1000
## semi parametric
ncc.np.sim <- replicate( n = Nsim, 
                         run_one_sim(sample.design = "ncc", 
                                     estimation.method = "NP"))

save(ncc.np.sim, file = "ncc.np.sim.Rdata")
