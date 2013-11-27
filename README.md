survMarkerTwoPhase
=============================================

This R package computes non-parametric and semi-parametric estimates of accuracy measures for risk prediction markers from survival data under two phase designs. For now, estimates can be obtained for the the case-cohort study design. Calculations for the nested case-control design will be implemented soon. 

Once sample weights are obtained, the function `survMTP.estimate` estimates the **AUC**, **TPR(c)**, **FPR(c)**, **PPV(c)**, and **NPV(c)** for for a specific timepoint and marker cutoff value c.

For more information, see references below. 


### Download and install the package



```r
# download the package from github
if (!require("devtools")) install.packages("devtools")
devtools::install_github("survMarkerTwoPhase", "mdbrown")
```



```r
library(survMarkerTwoPhase)
```

```
## Loading required package: survival
```

```
## Loading required package: splines
```

```r

# simulated data for illustration
data(SimData)
head(SimData)
```

```
##   survTime status        Y
## 1   0.1197      1  1.49310
## 2   1.0231      0 -0.73260
## 3   0.8282      0 -0.50211
## 4   2.0875      1  0.65758
## 5   4.6827      1  1.57806
## 6   0.3001      1  0.02419
```




### Case-Cohort Design

Estimation using a case-cohort subcohort design is permitted. Sample weights must first be calculated.

####Type I: 
 
Sample **n=150** from entire cohort, and include all participants with observed failures. 


```r
# generate a sub-cohort from SimData
set.seed(12321)
# create a sample index. 1 if sampled, 0 if not
N <- nrow(SimData)
sampleInd <- rep(0, N)

# sample all with observed failure time. (200 individuals)
sampleInd[SimData$status == 1] <- 1

# sample 150 more observations from the entire data set without
# replacement
sampleInd[sample(1:N, 150)] <- 1

table(sampleInd)  #total number of subcohort is 293 
```

```
## sampleInd
##   0   1 
## 207 293
```


Now we calculate the sample weights *w = 1/Pr(Sampled from cohort)* for each observation included in the sub-cohort. 

```r
# first calculate the Pr(Sampled from cohort) for each observation
sampleProb <- numeric(500)

# all non-censored observations were sampled, so their sample probability
# is 1
sampleProb[SimData$status == 1] <- 1
# all other individuals had a 150/N chance to be sampled
sampleProb[SimData$status == 0] <- 150/N

# the sample weights are 1/(probability of being sampled)
SimData$weights <- 1/sampleProb

subCohortData <- SimData[sampleInd == 1, ]
```

Here we estimate all measures available using non-parametric estimation methods. For non-parametric estimation, we need to provide the subcohort data consisting of status, survival time, marker, and weight, as well as survival time and status for the full cohort.   


```r
# estimate accuracy measures using non-parametric estimates by setting
# ESTmethod = 'NP' we have to supply subcohort data and full cohort data
survMTP.estimate(time = "survTime", event = "status", marker = "Y", weights = "weights", 
    subcohort.data = subCohortData, cohort.data = SimData, cohort.size = N, 
    estimation.method = "NP", study.design = "Case-Cohort", predict.time = 2, 
    cutpoint = 0)
```

```
## 
##  Non-parametric estimates under Case-Cohort study design.
## 
##         estimate
## AUC        0.749
## TPR(c)     0.699
## FPR(c)     0.394
## PPV(c)     0.331
## NPV(c)     0.878
## 
##  marker cutpoint: c = 0
```


Now estimate all measures  using semi-parametric estimation methods. These methods assume a Cox proportional hazards model. We only have to provide subcohort data in this case. 


```r
# estimate accuracy measures using semi-parametric estimates, we only need
# the sub-cohort data here
survMTP.estimate(time = "survTime", event = "status", marker = "Y", weights = "weights", 
    subcohort.data = subCohortData, cohort.size = N, estimation.method = "SP", 
    study.design = "Case-Cohort", predict.time = 2, cutpoint = 0)
```

```
## 
##  Semi-parametric estimates under Case-Cohort study design.
## 
##         estimate
## coef       1.088
## AUC        0.788
## TPR(c)     0.770
## FPR(c)     0.347
## PPV(c)     0.386
## NPV(c)     0.909
## 
##  marker cutpoint: c = 0
```


Only estimate the **AUC** and **TPR(0)**. 


```r
tmp <- survMTP.estimate(time = "survTime", event = "status", marker = "Y", weights = "weights", 
    subcohort.data = subCohortData, cohort.size = N, study.design = "Case-Cohort", 
    estimation.method = "SP", predict.time = 2, measures = c("AUC", "TPR"), 
    cutpoint = 0)
tmp
```

```
## 
##  Semi-parametric estimates under Case-Cohort study design.
## 
##         estimate
## coef       1.088
## AUC        0.788
## TPR(c)     0.770
## 
##  marker cutpoint: c = 0
```



```r
# access the estimates
tmp$estimates
```

```
##    coef    AUC    TPR
## 1 1.088 0.7876 0.7705
```


####Type II: 
 
In the second sampling design, we sample **n=100** from the strata of participants with observed failures and **n=100** from censored individuals. We need to account for this when we calculate the sampling weights. 


```r
# create a new sample index. 1 if sampled, 0 if not
sampleInd <- rep(0, N)

# sample 100 individuals with observed failure time.
CaseInd <- c(1:N)[SimData$status == 1]
sampleInd[sample(CaseInd, 100)] <- 1

# sample 100 individuals with censored failure time.
ControlInd <- c(1:N)[SimData$status == 0]
sampleInd[sample(ControlInd, 100)] <- 1

table(sampleInd)  #total number of subcohort is 200 
```

```
## sampleInd
##   0   1 
## 300 200
```


Calculate the sample weights.


```r
# first calculate the Pr(Sampled from cohort) for each observation
sampleProb <- numeric(500)

# all non-censored observations were sampled
sampleProb[SimData$status == 1] <- 100/length(CaseInd)
sampleProb[SimData$status == 0] <- 100/length(ControlInd)

SimData$weights2 <- 1/sampleProb

subCohortData2 <- SimData[sampleInd == 1, ]
```


Estimate all measures available, first using non-parametric estimation and then using semi-parametric estimates. 


```r
# Non-parametric estimates
survMTP.estimate(time = "survTime", event = "status", marker = "Y", weights = "weights2", 
    subcohort.data = subCohortData2, cohort.data = SimData, cohort.size = N, 
    estimation.method = "NP", study.design = "Case-Cohort", predict.time = 2, 
    cutpoint = 0)
```

```
## 
##  Non-parametric estimates under Case-Cohort study design.
## 
##         estimate
## AUC        0.694
## TPR(c)     0.590
## FPR(c)     0.377
## PPV(c)     0.340
## NPV(c)     0.822
## 
##  marker cutpoint: c = 0
```

```r

# Semi-parametric estimates
survMTP.estimate(time = "survTime", event = "status", marker = "Y", weights = "weights2", 
    subcohort.data = subCohortData2, cohort.size = N, estimation.method = "SP", 
    study.design = "Case-Cohort", predict.time = 2, cutpoint = 0)
```

```
## 
##  Semi-parametric estimates under Case-Cohort study design.
## 
##         estimate
## coef       0.908
## AUC        0.751
## TPR(c)     0.683
## FPR(c)     0.311
## PPV(c)     0.423
## NPV(c)     0.867
## 
##  marker cutpoint: c = 0
```

For more information see `?survMTP.estimate`. 


### References
Liu D, Cai T, Zheng Y. Evaluating the predictive value of biomarkers with stratified case-cohort design. *Biometrics* 2012, 4: 1219-1227.

Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














