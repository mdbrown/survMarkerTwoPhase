survMarkerTwoPhase
=============================================

This R package computes semi-parametric estimates of accuracy measures for risk prediction markers from survival data under two phase designs. For now, estimates can be obtained for the the case-cohort study design. Calculations for the nested case-control design will be implemented soon.

Once sample weights are obtained, the function `survMTP.estimate` estimates the $AUC$, $TPR(c)$, $FPR(c)$, $PPV(c)$, and $NPV(c)$ for for a specific timepoint and marker cutoff value c. Standard errors, and confidence intervals are also computed. 

For more information, see references below. 


### Download and install the package



```r
# download the package from github
if (!require("devtools")) install.packages("devtools")
```

```
## Loading required package: devtools
```

```r
devtools::install_github("survMarkerTwoPhase", "mdbrown")
```

```
## Installing github repo(s) survMarkerTwoPhase/master from mdbrown
## Downloading survMarkerTwoPhase.zip from
## https://github.com/mdbrown/survMarkerTwoPhase/archive/master.zip
## Installing package from
## C:\Users\mdbrown\AppData\Local\Temp\RtmpgT33AL/survMarkerTwoPhase.zip
## Installing survMarkerTwoPhase "C:/PROGRA~1/R/R-30~1.1/bin/x64/R" --vanilla
## CMD INSTALL \
## "C:\Users\mdbrown\AppData\Local\Temp\RtmpgT33AL\survMarkerTwoPhase-master"
## \ --library="C:/Users/mdbrown/Documents/R/win-library/3.0" \
## --with-keep.source --install-tests
```

```r


library(survMarkerTwoPhase)
```

```
## Loading required package: survival Loading required package: splines
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

Estimation using a case-cohort subcohort design is permitted. Sample weights must first be calculated based.

####Type I: 
 
Sample $n=150$ from entire cohort, and include all participants with observed failures. 


```r
# generate a sub-cohort from SimData
set.seed(12321)
# create a sample index. 1 if sampled, 0 if not
N <- nrow(SimData)
sampleInd <- rep(0, N)

# sample all with observed failure time. (200 individuals)
sampleInd[SimData$status == 1] <- 1

# sample 150 more observations from the entire data set without replacement
sampleInd[sample(1:N, 150)] <- 1

table(sampleInd)  #total number of subcohort is 293 
```

```
## sampleInd
##   0   1 
## 207 293
```


Now we calculate the sample weights $w_i = 1/Pr(Sampled from cohort)$ for each observation included in the sub-cohort. 

```r

# first calculate the Pr(Sampled from cohort) for each observation
sampleProb <- numeric(500)

# all non-censored observations were sampled, so their sample probability is
# 1
sampleProb[SimData$status == 1] <- 1
sampleProb[SimData$status == 0] <- 150/N

SimData$weights <- 1/sampleProb

subCohortData <- SimData[sampleInd == 1, ]
```


Here we estimate all measures available. 


```r
# estimate accuracy measures using only the subcohort data
survMTP.estimate(time = subCohortData$survTime, event = subCohortData$status, 
    marker = subCohortData$Y, weights = subCohortData$weights, study.design = "Case-Cohort", 
    predict.time = 2, cutpoint = 0)
```

```
## 
## Semi-parametric estimates under Case-Cohort study design.
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.088     0.167         0.761       1.414 
## AUC        0.788     0.034         0.713       0.847 
## TPR(c)     0.770     0.058         0.639       0.864 
## FPR(c)     0.347     0.057         0.245       0.465 
## PPV(c)     0.386     0.042         0.308       0.471 
## NPV(c)     0.909     0.020         0.862       0.942 
## 
##  marker cutpoint: c = 0
```


Only estimate the $AUC$ and $TPR(0)$. 


```r
tmp <- survMTP.estimate(time = subCohortData$survTime, event = subCohortData$status, 
    marker = subCohortData$Y, weights = subCohortData$weights, study.design = "Case-Cohort", 
    predict.time = 2, measures = c("AUC", "TPR"), cutpoint = 0)
tmp
```

```
## 
## Semi-parametric estimates under Case-Cohort study design.
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.088     0.167         0.761       1.414 
## AUC        0.788     0.034         0.713       0.847 
## TPR(c)     0.770     0.058         0.639       0.864 
## 
##  marker cutpoint: c = 0
```



```r
# access the estimates
tmp$estimates
```

```
##   coef    AUC    TPR 
## 1.0876 0.7876 0.7705
```

```r

# and the confidence bounds
tmp$CIbounds
```

```
##         coef   AUC    TPR
## upper 1.4140 0.847 0.8644
## lower 0.7613 0.713 0.6386
```


####Type II: 
 
In the second sampling design, we sample $n=100$ from the strata of participants with observed failures and $n=100$ from censored individuals. We need to account for this when we calculate the sampling weights. 


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


Estimate all measures available. 


```r
# estimate accuracy measures using only the subcohort data
survMTP.estimate(time = subCohortData2$survTime, event = subCohortData2$status, 
    marker = subCohortData2$Y, weights = subCohortData2$weights2, study.design = "Case-Cohort", 
    predict.time = 2, cutpoint = 0)
```

```
## 
## Semi-parametric estimates under Case-Cohort study design.
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       0.908     0.135         0.644       1.171 
## AUC        0.751     0.041         0.661       0.823 
## TPR(c)     0.683     0.090         0.489       0.829 
## FPR(c)     0.311     0.083         0.175       0.490 
## PPV(c)     0.423     0.060         0.312       0.544 
## NPV(c)     0.867     0.028         0.802       0.912 
## 
##  marker cutpoint: c = 0
```







For more information see `?survMTP.estimate`. 



### References
Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














