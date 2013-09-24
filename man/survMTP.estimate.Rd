\name{survMTP.estimate}
\alias{survMTP.estimate}
\title{Semi-parametric estimates of accuracy measures for a risk prediction marker from survival data under two phase designs}
\usage{
  survEstMeasures(time, event, marker, weights = NULL, study.design = "Case-Cohort", 
    predict.time,
    measures = c("all"), cutpoint = median(marker),
    CImethod = "logit.transformed", alpha = 0.05)
}
\arguments{
  \item{time}{numeric vector, time to event variable}

  \item{event}{numeric vector, indicator for the status of
  event of interest. event = 0 for censored observations,
  and event = 1 for event of interest.}

  \item{marker}{numeric vector, marker variable of
  interest}
  
  \item{weights}{numeric vector of sample weights from a case-cohort sample design. Should be equal to 1/Pr(sampled from cohort). }
    \item{study.design}{ String indicating the type of study design. Only "Case-Cohort" is available now, "Nested Case-Control" will be added soon. }
  \item{predict.time}{numeric value of the timepoint of
  interest for which to estimate the risk measures}

  \item{measures}{character vector that can contain the
  following 'AUC', 'TPR', 'FPR', 'PPV', 'NPV' or 'all'
  (default).}

  \item{cutpoint}{numeric value indicating the value of the
  cutpoint 'c' at which to estimate 'FPR', 'TPR', 'NPV' or
  'PPV'. default is \code{median(marker)}.}

  \item{CImethod}{character string of either
  'logit.transformed' (default) or 'standard' indicating
  whether normal approximated confidence intervals should
  be calculated using logistic transformed values or the
  standard method.}

  \item{alpha}{alpha value for confidence intervals.
  (1-alpha)*100% confidence intervals are provided. default
  is alpha = 0.05.}
}
\value{
  a list with components \item{estimates}{point estimates
  for risk measures} \item{se}{standard errors for
  estimates} \item{CIbounds}{bounds for (1-alpha)*100%
  confidence interval} \item{model.fit}{object of type
  'coxph' from fitting the model \code{coxph(Surv(time,
  event)~Y)}} \item{cutoff, CImethod, SEmethod,
  predict.time, alpha}{ function inputs}
}
\description{
  This function estimates the AUC, TPR(c), FPR(c), PPV(c),
  and NPV(c) for for a specific timepoint and marker cutoff
  value c. Standard errors, and confidence intervals are
  also computed. Either analytic or bootstrap standard
  errors can be computed.
}
\examples{
data(SimData)

tmp <- survEstMeasures(time =SimData$survTime, event = SimData$status, marker = SimData$Y, predict.time = 2, cutpoint = 0)
tmp

survEstMeasures(time =SimData$survTime, event = SimData$status, marker = SimData$Y, predict.time = 2, measures = c("AUC", "TPR"), SEmethod = 'bootstrap', bootstraps = 50, cutpoint = 0)

#generate a sub-cohort from SimData
set.seed(12321)
#create a sample index. 1 if sampled, 0 if not
N <- nrow(SimData)
sampleInd <- rep(0, N)

# sample all with observed failure time. (200 individuals)
sampleInd[SimData$status==1] <- 1

#sample 150 more observations from the entire data set without replacement
sampleInd[sample(1:N, 150)] <- 1

table(sampleInd) #total number of subcohort is 293 

## now calculate sample weights
#first calculate the Pr(Sampled from cohort) for each observation
sampleProb <- numeric(500)
#all non-censored observations were sampled, so their sample probability is 1
sampleProb[SimData$status==1] <- 1 
sampleProb[SimData$status==0] <- 150/N

SimData$weights <- 1/sampleProb

subCohortData <- SimData[sampleInd==1,]

#estimate accuracy measures using only the subcohort data
survEstMeasures(time =subCohortData$survTime, 
                event = subCohortData$status, 
                marker = subCohortData$Y,
                weights = subCohortData$weights, 
                predict.time = 2, 
                 cutpoint = 0, SEmethod = "normal")



}
\references{
  Liu D, Cai T, Zheng Y. Evaluating the predictive value of biomarkers
  with stratified case-cohort design. *Biometrics* 2012, 4: 1219-1227.

  Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance
  of markers for future events. Lifetime Data Analysis.
  2008, 14: 86-113.

  Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent
  predictive values of prognostic biomarkers with failure
  time outcome. JASA 2008, 103: 362-368.
}

