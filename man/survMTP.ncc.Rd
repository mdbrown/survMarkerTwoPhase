\name{survMTP.ncc}
\alias{survMTP.ncc}
\title{Accuracy measures for a risk prediction marker under nested case-control study design}
\usage{
  survMTP.ncc(time, event, marker, 
                           subcoh, 
                           id, 
                           data,
                           risk.sets, 
                           estimation.method = "NP", 
                           predict.time,
                           alpha = 0.05, 
                           Npert = 500, 
                           marker.cutpoint = 'median' )
}
\arguments{
  \item{time}{time to event}

  \item{event}{indicator for the status of
  event of interest. event = 0 for censored observations,
  and event = 1 for event of interest.}

  \item{marker}{marker variable of
  interest}
  
    \item{subcoh}{indicator for subjects included in the subcohort sample (1=included, 0=not included)}
    
    \item{id}{ unique numeric identifier for each observation.}
  
  \item{data}{data frame in which to look for input variables. This should be the full cohort data, including time, event, subcoh, id, and marker values when subcoh = 1. }
  
  \item{risk.sets}{ a data frame or matrix identifying risk sets. Specifically this should have one row for each case (event==1), consisting of the case id, followed by corresponding matched id's. See below for a more concrete example. 
  }
  
      \item{estimation.method}{should non-parametric ('NP') or semi-parametric ('SP') estimates be calculated. Semi-parametric methods assume a proportional hazards model.}
  
  \item{predict.time}{numeric value of the timepoint of
  interest for which to estimate the risk measures}

  \item{alpha}{alpha value for confidence intervals.
  (1-alpha)*100\% confidence intervals are provided. default
  is alpha = 0.05.}

 \item{Npert}{ number of perturbations to use when calculating standard errors. Default is 500.}

  \item{marker.cutpoint}{numeric value indicating the value of the
  cutpoint 'c' at which to estimate 'FPR', 'TPR', 'NPV' or
  'PPV'. default is \code{'median'} which takes cutpoint as the marker median.}

}
\value{
  a list with components 
  \item{estimates}{point estimates
  for risk measures}
    \item{se}{standard errors for 
   risk measure estimates} 
       \item{CIbounds}{ (1-alpha)*100 \% confidence bounds for measures.
       } 
  \item{model.fit}{*only returned if estimation.method = 'SP'* object of type
  'coxph' from fitting the model 
    \code{coxph(Surv(time, event)~Y)}}
  \item{marker.cutpoint, estimation.method, predict.time}{ function inputs}
}
\description{
  This function estimates the AUC, TPR(c), FPR(c), PPV(c),
  and NPV(c) for for a specific timepoint and marker cutpoint
  value c with survival data from a nested case-control subcohort design. 
}
\examples{
#simulated data for illustration 
data(SimData)

# use the function 'ccwc' from the 
# epi package to generate a subcohort. 
require('Epi')

#2 matched controls for each case
nmatch = 2

cohortData_ncc <- SimData
cohortData_ncc$id = 1:dim(cohortData_ncc)[1]

subcohort_ncc <- ccwc(exit = survTime, fail= status, data = cohortData_ncc, controls = nmatch) # match 2 controls for each case. 

#indicator for inclusion in the subcohort
sampleInd = rep(0, nrow(cohortData_ncc)) #initialize all equal to zero
sampleInd[subcohort_ncc$Map] = 1

cohortData_ncc$subcohort = sampleInd  

table(sampleInd) #subcohort sample size of 352

#marker data is unavailable for those not in the subcohort
cohortData_ncc$Y[sampleInd==0] = NA  

#now we need to build the risk set matrix,
#which will be dimension (# of cases) x (nmatch + 1), so 200x3 here
#each row denotes a risk set, with the case id followed by the matched control ids
RiskSets <- matrix(nrow = sum(cohortData_ncc$status), ncol = nmatch +1)

for(i in subcohort_ncc$Set){
  RiskSets[i, ] <- unlist(subset(subcohort_ncc, Set==i, select=Map))
}

#Nonparametric estimates 
survMTP.ncc(time = survTime, 
         event = status, 
         marker=Y, 
         subcoh=subcohort, 
         id = id, 
         data = cohortData_ncc,
         risk.sets = RiskSets, 
         estimation.method = "NP", 
         predict.time = 2 , 
         marker.cutpoint = 0)

#Semiparametric estimates 
survMTP.ncc(time = survTime, 
         event = status, 
         marker=Y, 
         subcoh=subcohort, 
         id = id, 
         data = cohortData_ncc,
         risk.sets = RiskSets, 
         estimation.method = "SP", 
         predict.time = 2 , 
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

