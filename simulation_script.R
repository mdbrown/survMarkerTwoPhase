## Script to check the standard error estimates under CCH designs

#need the function to simulate data
source("../../power/survAccuracyMeasuresPower/subroutines.R")

#get betas needed for AUC of .5, .6, .7, .8, .9

round(get.Betas(parameter = "AUC", a = .1, predict.time =2, cutoff = 0, 
                parval.H0 =.7, parval.Ha=.9,  f_x = dnorm ), 4)
  
AUCVec  = c(.5, .6, .7, .8, .9)
betaVec = c(0, .3247, 0.6759, 1.1220, 1.9387)

NSim = 1000
N <- 500
#type I and type II estimates 
type1.estimates <- data.frame("AUC" = numeric(NSim*length(betaVec)),
                              "AUCest" = numeric(NSim*length(betaVec)), 
                              "AUCSE" = numeric(NSim*length(betaVec)),
                              "beta" = numeric(NSim*length(betaVec)),
                              "betaest"= numeric(NSim*length(betaVec)), 
                              "betaSE" = numeric(NSim*length(betaVec)), 
                              "subN" = numeric(NSim*length(betaVec)))
type2.estimates <- data.frame("AUC" = numeric(NSim*length(betaVec)),
                              "AUCest" = numeric(NSim*length(betaVec)), 
                              "AUCSE" = numeric(NSim*length(betaVec)),
                              "beta" = numeric(NSim*length(betaVec)),
                              "betaest"= numeric(NSim*length(betaVec)), 
                              "betaSE" = numeric(NSim*length(betaVec)))
cohort.estimates <- data.frame("AUC" = numeric(NSim*length(betaVec)),
                              "AUCest" = numeric(NSim*length(betaVec)), 
                              "AUCSE" = numeric(NSim*length(betaVec)),
                              "beta" = numeric(NSim*length(betaVec)),
                              "betaest"= numeric(NSim*length(betaVec)), 
                              "betaSE" = numeric(NSim*length(betaVec)))

index = 0
for(j in 1:length(AUCVec)){

for( i in 1:NSim){
  index = index + 1; 
#simulate cohorts with 40% censoring, 500 individuals
SimData <- SIM.data.singleMarker(nn = N, 
           beta = betaVec[j], 
           cens.perc = .7)
  
  tmp <- survAM.estimate(time =SimData$xi, 
                          event = SimData$di, 
                          marker = SimData$Y,
                          measures = "AUC",
                          predict.time = 2, 
                          cutpoint = 0)
 cohort.estimates$AUC[index] = AUCVec[j]
 cohort.estimates$AUCest[index] = tmp$estimates[2]
 cohort.estimates$SE[index] = tmp$se[2]
  
  cohort.estimates$beta[index] = betaVec[j]
  cohort.estimates$betaest[index] <- tmp$estimates[1]
  cohort.estimates$betaSE[index] <- tmp$se[1]
#CCH type I  

sampleInd <- rep(0, N)

# sample all with observed failure time. (200 individuals)
sampleInd[SimData$di==1] <- 1

#sample 150 more observations from the entire data set without replacement
sampleInd[sample(1:N, 100)] <- 1
sampleProb <- numeric(N)

#all non-censored observations were sampled, so their sample probability is 1
sampleProb[SimData$di==1] <- 1 
sampleProb[SimData$di==0] <- 100/N

SimData$weights <- 1/sampleProb

subdat.type1 <- SimData[sampleInd==1,]

#estimate accuracy measures using only the subcohort data
tmp <- survMTP.estimate(time =subdat.type1$xi, 
                 event = subdat.type1$di, 
                 marker = subdat.type1$Y,
                 weights = subdat.type1$weights,
                 study.design = "Case-Cohort",
                 measures = "AUC",
                 predict.time = 2, 
                 cutpoint = 0)

type1.estimates$AUC[index] = AUCVec[j]
type1.estimates$AUCest[index] = tmp$estimates[2]
type1.estimates$SE[index] = tmp$se[2]
type1.estimates$subN[index] = sum(sampleInd)
  
  type1.estimates$beta[index] = betaVec[j]
  type1.estimates$betaest[index] <- tmp$estimates[1]
  type1.estimates$betaSE[index] <- tmp$se[1]
#### Type 2

sampleInd <- rep(0, N)

# sample 100 individuals with observed failure time. 
CaseInd <- c(1:N)[SimData$di==1] 
sampleInd[sample(CaseInd, 100)] <- 1

# sample 100 individuals with censored failure time. 
ControlInd <- c(1:N)[SimData$di==0] 
sampleInd[sample(ControlInd, 100)] <- 1

sampleProb <- numeric(500)

#all non-censored observations were sampled
sampleProb[SimData$di==1] <- 100/length(CaseInd) 
sampleProb[SimData$di==0] <- 100/length(ControlInd)

SimData$weights2 <- 1/sampleProb

subdat.type2 <- SimData[sampleInd==1,]

#estimate accuracy measures using only the subcohort data
tmp <- survMTP.estimate(time =subdat.type2$xi, 
                        event = subdat.type2$di, 
                        marker = subdat.type2$Y,
                        weights = subdat.type2$weights2,
                        study.design = "Case-Cohort",
                        measures = "AUC",
                        predict.time = 2, 
                        cutpoint = 0)

type2.estimates$AUC[index] = AUCVec[j]
type2.estimates$AUCest[index] = tmp$estimates[2]
type2.estimates$SE[index] = tmp$se[2]
  type2.estimates$beta[index] = betaVec[j]
  type2.estimates$betaest[index] <- tmp$estimates[1]
  type2.estimates$betaSE[index] <- tmp$se[1]
print(paste(AUCVec[j], i, sep = ":"))

}
}

ddply(type1.estimates, "AUC", function(df) c(mean(df$AUCest), mean(df$SE), sd(df$AUCest), 
                            mean(df$beta),  mean(df$betaest), mean(df$betaSE), sd(df$betaest)))
ddply(type2.estimates, "AUC", function(df) c(mean(df$AUCest), mean(df$SE), sd(df$AUCest)))
ddply(cohort.estimates, "AUC", function(df)c(mean(df$AUCest), mean(df$SE), sd(df$AUCest), 
                                             mean(df$beta),  mean(df$betaest), mean(df$betaSE), sd(df$betaest)))





