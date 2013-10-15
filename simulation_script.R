## Script to check the standard error estimates under CCH designs

#need the function to simulate data
source("../../power/survAccuracyMeasuresPower/subroutines.R")
source("../PowerSAM/subroutines.R")
#get betas needed for AUC of .5, .6, .7, .8, .9

round(get.Betas(parameter = "AUC", a = .1, predict.time =2, cutoff = 0, 
                parval.H0 =.75, parval.Ha=.9,  f_x = dnorm ), 4)
  
AUCVec  = c(.5, .6, .7, .8, .9)
betaVec = c(0, .3247, 0.6759, 1.1220, 1.9387)


mybeta <- 0.8789

NSim = 1000
N <- 1000
#type I and type II estimates 
cch.1.est <- matrix(ncol = 6, nrow = NSim)
cch.1.se  <- matrix(ncol = 6, nrow = NSim)
cch.1.N <- numeric(NSim)

cch.2.est <- matrix(ncol = 6, nrow = NSim)
cch.2.se  <- matrix(ncol = 6, nrow = NSim)

cohort.est <- matrix(ncol = 6, nrow = NSim)
cohort.se  <- matrix(ncol = 6, nrow = NSim)


setwd(12321)
index = 0

for( i in 1:NSim){
  index = index + 1; 
#simulate cohorts with 40% censoring, 500 individuals
SimData <- SIM.data.singleMarker(nn = N, 
           beta = mybeta, 
           cens.perc = .7)
  
  tmp <- survAM.estimate(time =SimData$xi, 
                          event = SimData$di, 
                          marker = SimData$Y,
                          predict.time = 2, 
                          cutpoint = 0)
 cohort.est[i,] = tmp$estimates
 cohort.se[i,]  = tmp$se

#CCH type I  

sampleInd <- rep(0, N)

# sample all with observed failure time. (300 individuals)
sampleInd[SimData$di==1] <- 1

#sample 150 more observations from the entire data set without replacement
sampleInd[sample(1:N, 300)] <- 1
sampleProb <- numeric(N)

#all non-censored observations were sampled, so their sample probability is 1
sampleProb[SimData$di==1] <- 1 
sampleProb[SimData$di==0] <- 300/N

SimData$weights <- 1/sampleProb

subdat.type1 <- SimData[sampleInd==1,]

#estimate accuracy measures using only the subcohort data
tmp <- survMTP.estimate(time =subdat.type1$xi, 
                 event = subdat.type1$di, 
                 marker = subdat.type1$Y,
                 weights = subdat.type1$weights,
                        cohortN = N, 
                 study.design = "Case-Cohort",
                 predict.time = 2, 
                 cutpoint = 0)

  cch.1.est[i,] = tmp$estimates
  cch.1.se[i,]  = tmp$se
  cch.1.N[i] = sum(sampleInd)
  #### Type 2

sampleInd <- rep(0, N)

# sample 100 individuals with observed failure time. 
CaseInd <- c(1:N)[SimData$di==1] 
sampleInd[sample(CaseInd, 250)] <- 1

# sample 100 individuals with censored failure time. 
ControlInd <- c(1:N)[SimData$di==0] 
sampleInd[sample(ControlInd, 250)] <- 1

sampleProb <- numeric(500)

#all non-censored observations were sampled
sampleProb[SimData$di==1] <- 250/length(CaseInd) 
sampleProb[SimData$di==0] <- 250/length(ControlInd)

SimData$weights2 <- 1/sampleProb

subdat.type2 <- SimData[sampleInd==1,]

#estimate accuracy measures using only the subcohort data
tmp <- survMTP.estimate(time =subdat.type2$xi, 
                        event = subdat.type2$di, 
                        marker = subdat.type2$Y,
                        weights = subdat.type2$weights2,
                        cohortN = N, 
                        study.design = "Case-Cohort",
                        predict.time = 2, 
                        cutpoint = 0)

  cch.2.est[i,] = tmp$estimates
  cch.2.se[i,]  = tmp$se

print(paste( i))

}
cohort.est <- as.data.frame(cohort.est)
cohort.se <- as.data.frame(cohort.se)

cch.1.est <- as.data.frame(cch.1.est)
cch.1.se <- as.data.frame(cch.1.se)

cch.2.est <- as.data.frame(cch.2.est)
cch.2.se <- as.data.frame(cch.2.se)

names(cch.1.est) <- names(cch.1.se) <- names(cch.2.est) <- names(cch.2.se) <- names(cohort.se) <- names(cohort.est) <- names(tmp$estimates)



## cohort
colMeans(cohort.est)
round(apply(cohort.est, 2, sd) -colMeans(cohort.se), 4)

#type 1
colMeans(cch.1.est)
round(apply(cch.1.est, 2, sd) -colMeans(cch.1.se), 4)

#type 1
colMeans(cch.2.est, na.rm = TRUE)
round(apply(cch.2.est, 2, sd) -colMeans(cch.2.se), 4)

#make a little table for each 















load("simOutput.Rdata")

ddply(type1.estimates, "AUC", function(df)c("MeanAUCest" = mean(df$AUCest), 
                                            "MeanAUCSE" = mean(df$SE), 
                                            "EmpAUCSE"  = sd(df$AUCest), 
                                            "beta" = mean(df$beta),  
                                            "MeanbetaEST" = mean(df$betaest), 
                                            "MeanbetaSE" = mean(df$betaSE), 
                                            "EmpbetaSE" = sd(df$betaest)))

AUC MeanAUCest   MeanAUCSE    EmpAUCSE   beta  MeanbetaEST MeanbetaSE  EmpbetaSE
1 0.5  0.4996356 0.024684076 0.025140983 0.0000 0.0009909321 0.07892750 0.08205906
2 0.6  0.5992227 0.023699374 0.023975690 0.3247 0.3294330829 0.07993058 0.08161297
3 0.7  0.6970220 0.020989044 0.022175399 0.6759 0.6780287405 0.08128971 0.08745908
4 0.8  0.7974996 0.015663825 0.016393401 1.1220 1.1319576301 0.08445859 0.09303013
5 0.9  0.8954575 0.008856503 0.008975003 1.9387 1.9497646535 0.09875168 0.10606743

ddply(type2.estimates, "AUC", function(df)c("MeanAUCest" = mean(df$AUCest), 
                                            "MeanAUCSE" = mean(df$SE), 
                                            "EmpAUCSE"  = sd(df$AUCest), 
                                            "beta" = mean(df$beta),  
                                            "MeanbetaEST" = mean(df$betaest), 
                                            "MeanbetaSE" = mean(df$betaSE), 
                                            "EmpbetaSE" = sd(df$betaest)))
AUC MeanAUCest   MeanAUCSE    EmpAUCSE   beta MeanbetaEST MeanbetaSE  EmpbetaSE
1 0.5  0.5000122 0.024210693 0.024425194 0.0000 0.002818663 0.07756309 0.07975210
2 0.6  0.5977038 0.023389232 0.023178463 0.3247 0.324607724 0.07875714 0.07821797
3 0.7  0.6976133 0.020921332 0.021488859 0.6759 0.679432705 0.08072867 0.08205567
4 0.8  0.7981976 0.015923229 0.016440535 1.1220 1.129830649 0.08518855 0.09043359
5 0.9  0.8962485 0.009183705 0.009503245 1.9387 1.946625503 0.10322186 0.10839158

ddply(cohort.estimates, "AUC", function(df)c("MeanAUCest" = mean(df$AUCest), 
                                             "MeanAUCSE" = mean(df$SE), 
                                             "EmpAUCSE"  = sd(df$AUCest), 
                                             "beta" = mean(df$beta),  
                                             "MeanbetaEST" = mean(df$betaest), 
                                             "MeanbetaSE" = mean(df$betaSE), 
                                             "EmpbetaSE" = sd(df$betaest)))


AUC MeanAUCest   MeanAUCSE   EmpAUCSE   beta  MeanbetaEST MeanbetaSE  EmpbetaSE
1 0.5  0.4995863 0.018013997 0.01797672 0.0000 0.0003545359 0.05774681 0.05793222
2 0.6  0.5989182 0.017564602 0.01731011 0.3247 0.3241587557 0.05849438 0.05769291
3 0.7  0.6980886 0.016121246 0.01699209 0.6759 0.6735228188 0.06102150 0.06379761
4 0.8  0.7995466 0.012828753 0.01342071 1.1220 1.1263176713 0.06791604 0.07256237
5 0.9  0.8984744 0.007675964 0.00784358 1.9387 1.9428456342 0.08776859 0.09003101


