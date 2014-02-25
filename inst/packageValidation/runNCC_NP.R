library(survMarkerTwoPhase)

for(arg in commandArgs(TRUE)) eval(parse(text=arg)) #expecting batch


new_counter <- function() {
  i <- 0
  function() {
    i <<- i + 1
    i
  }
}


source("simulateData.R")
Nsim = 100
## semi parametric
mycounter <- new_counter()

ncc.np.sim <- replicate( n = Nsim, 
                         run_one_sim(sample.design = "ncc", 
                                     estimation.method = "NP", 
                                     counter= mycounter))

save(ncc.np.sim, file = paste("ncc.np.sim_", batch, ".Rdata", sep=""))
