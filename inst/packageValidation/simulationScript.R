# script to run all the simulations. calls the functions from 'simulateData.R' in this directory
source("simulateData.R")

Nsim = 1000

######
## cch 
##
## non parametric
set.seed(1221)
cch.np.sim <- replicate( n = Nsim, 
                         run_one_sim(sample.design = "cch", 
                                     estimation.method = "NP"))

save(cch.np.sim, file = "cch.np.sim.Rdata")

## semi parametric
cch.sp.sim <- replicate( n = Nsim, 
                         run_one_sim(sample.design = "cch", 
                                     estimation.method = "SP"))

save(cch.sp.sim, file = "cch.sp.sim.Rdata")
                                                  


######
## ncc
##
## non parametric
ncc.np.sim <- replicate( n = Nsim, 
                         run_one_sim(sample.design = "ncc", 
                                     estimation.method = "NP"))

save(ncc.np.sim, file = "ncc.np.sim.Rdata")

## semi parametric
ncc.sp.sim <- replicate( n = Nsim, 
                         run_one_sim(sample.design = "ncc", 
                                     estimation.method = "SP"))

save(ncc.sp.sim, file = "ncc.sp.sim.Rdata")
