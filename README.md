survMarkerTwoPhase
=============================================

This R package computes non-parametric and semi-parametric estimates of accuracy measures for risk prediction markers from survival data under two phase study designs, namely the case-cohort and the nested case-control study designs.

The accuracy measures that can be estimated include the **AUC**, **TPR(c)**, **FPR(c)**, **PPV(c)**, and **NPV(c)** for for a specific timepoint and fixed marker cutoff value c.

Below is a brief tutorial to get you started. For more information regarding estimation procedures see the references below. 


#Tutorial

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

For the case-cohort design we sample **n=150** from entire cohort, and include all participants with observed failures. 


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

To estimate accuracy measures for the case-cohort design, sample weights must be calculated. The sample weights *w = 1/Pr(Sampled from cohort)* for each observation included in the sub-cohort. 

```r
# first calculate the Pr(Sampled from cohort) for each observation
sampleProb <- numeric(500)

# all non-censored observations were sampled, so their sample probability is
# 1
sampleProb[SimData$status == 1] <- 1
# all other individuals had a 150/N chance to be sampled
sampleProb[SimData$status == 0] <- 150/N

# the sample weights are 1/(probability of being sampled)
SimData$weights <- 1/sampleProb

# indicator of inclusion in the subcohort
SimData$subcohort = sampleInd
```

Here we estimate using non-parametric estimation methods.   


```r
# estimate accuracy measures using non-parametric estimates by setting
# ESTmethod = 'NP'
survMTP.cch(time = survTime, event = status, marker = Y, weights = weights, 
    subcoh = subcohort, data = SimData, estimation.method = "NP", predict.time = 2, 
    marker.cutpoint = 0)
```

```
## Called from: survMTP.cch(time = survTime, event = status, marker = Y, weights = weights, 
##     subcoh = subcohort, data = SimData, estimation.method = "NP", 
##     predict.time = 2, marker.cutpoint = 0)
```

```
## Warning: the condition has length > 1 and only the first element will be
## used
```

```
## Error: row names must be 'character' or 'integer', not 'double'
```


Now estimate measures  using semi-parametric estimation methods. These methods assume a Cox proportional hazards model.


```r
# estimate accuracy measures using semi-parametric estimates, we only need
# the sub-cohort data here
survMTP.cch(time = survTime, event = status, marker = Y, weights = weights, 
    subcoh = subcohort, data = SimData, estimation.method = "SP", predict.time = 2, 
    marker.cutpoint = 0)
```

```
## Called from: survMTP.cch(time = survTime, event = status, marker = Y, weights = weights, 
##     subcoh = subcohort, data = SimData, estimation.method = "SP", 
##     predict.time = 2, marker.cutpoint = 0)
```

```
## Error: object 'cutpoint' not found
```



For more information see `?survMTP.estimate`. 





### References
Liu D, Cai T, Zheng Y. Evaluating the predictive value of biomarkers with stratified case-cohort design. *Biometrics* 2012, 4: 1219-1227.

Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














