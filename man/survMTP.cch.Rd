\name{survMTP.cch}
\alias{survMTP.cch}
\title{Accuracy measures for a risk prediction marker under case-cohort study design}
\usage{
  survMTP.cch(time, event, marker, 
                             weights, 
                             subcoh,
                             data, 
                             estimation.method = 'NP', 
                             predict.time, 
                             marker.cutpoint = 'median')
}
\arguments{
  \item{time}{time to event}

  \item{event}{indicator for the status of
  event of interest. event = 0 for censored observations,
  and event = 1 for event of interest.}

  \item{marker}{marker variable of
  interest}
  
  \item{weights}{sample weights from a case-cohort sample design. Should be equal to 1/Pr(sampled from cohort). }
  
    \item{subcoh}{indicator for subjects included in the subcohort sample (1=included, 0=not included)}
  
  \item{data}{data frame in which to look for input variables. This should be the full cohort data, including time, event, weights, subcoh, and marker values when subcoh = 1. }
  
      \item{estimation.method}{should non-parametric ('NP') or semi-parametric ('SP') estimates be calculated. Semi-parametric methods assume a proportional hazards model.}
  
  \item{predict.time}{numeric value of the timepoint of
  interest for which to estimate the risk measures}

  \item{marker.cutpoint}{numeric value indicating the value of the
  cutpoint 'c' at which to estimate 'FPR', 'TPR', 'NPV' or
  'PPV'. default is \code{'median'} which takes cutpoint as the marker median.}

}
\value{
  a list with components 
  \item{estimates}{point estimates
  for risk measures} 
  \item{model.fit}{*only returned if estimation.method = 'SP'* object of type
  'coxph' from fitting the model 
    \code{coxph(Surv(time, event)~Y)}}
  \item{marker.cutpoint, estimation.method, predict.time}{ function inputs}
}
\description{
  This function estimates the AUC, TPR(c), FPR(c), PPV(c),
  and NPV(c) for for a specific timepoint and marker cutoff
  value c with survival data from a case-cohort subcohort design. 
}
\examples{

#simulated data for illustration
data(SimData)

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

cohortData_cch <- SimData

#first calculate the Pr(Sampled from cohort) for each observation
sampleProb <- numeric(500)

#all non-censored observations were sampled, so their sample probability is 1
sampleProb[cohortData_cch$status==1] <- 1 
#all other individuals had a 150/N chance to be sampled 
sampleProb[cohortData_cch$status==0] <- 150/N

#the sample weights are 1/(probability of being sampled)
cohortData_cch$weights <- 1/sampleProb

#indicator of inclusion in the subcohort
cohortData_cch$subcohort = sampleInd   

#estimate accuracy measures using non-parametric estimates
#by setting estimation.method = "NP"
survMTP.cch(time =survTime, 
            event = status, 
            marker = Y,
            weights = weights,
            subcoh = subcohort,
            data = cohortData_cch, 
            estimation.method = "NP",  
            predict.time = 2, 
            marker.cutpoint = 0)

#estimate accuracy measures using semi-parametric estimates
survMTP.cch(time =survTime, 
            event = status, 
            marker = Y,
            weights = weights,
            subcoh = subcohort,
            data = cohortData_cch, 
            estimation.method = "SP",  
            predict.time = 2, 
            marker.cutpoint = 0)


}
\references{
1. Liu D, Cai T, Zheng Y. Evaluating the predictive value of biomarkers with stratified case-cohort design. *Biometrics* 2012, 4: 1219-1227.

2. Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

3. Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.

4. Cai T. and Zheng Y . Resampling Procedures for Making Inference under Nested Case-control Studies. *JASA* 2013 (in press). 

5. Cai T and Zheng Y/*, Evaluating Prognostic Accuracy of Biomarkers under nested case-control studies. *Biostatistics* 2012, 13,1, 89-100. 

(* equal contributor and corresponding author). 

}

