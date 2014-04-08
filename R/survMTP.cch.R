 
survMTP.cch <- function(time, event, marker, weights, 
                             subcoh,
                             data, 
                             estimation.method = 'NP', 
                             predict.time, 
                             marker.cutpoint = 'median',
                             ci.method = "logit.transformed",
                             alpha=0.05
                           ){
  
  
  # checks
  stopifnot(is.data.frame(data))
  
  time <- eval(substitute(time), data)
  event <- 1*eval(substitute(event), data)
  marker <- eval(substitute(marker), data)
  weights <- eval(substitute(weights), data)
  vi <- 1*eval(substitute(subcoh), data)

  stopifnot(is.element(estimation.method, c("NP", "SP")))
  stopifnot(is.numeric(predict.time))
  if(marker.cutpoint=='median') marker.cutpoint =  median(eval(substitute(marker), data))
  stopifnot(is.numeric(marker.cutpoint))
  
  #set some defaults
  measures = c('all')
  cutoff <- marker.cutpoint
  cutoff.type = "none"; #cutoffN = 100 #functionality to be added later
  SEmethod ="normal"; #bootstraps = 10

  

  if(length(cutoff)==0) cutoff = NA;  
  
  #data.frame checks
  if(!is.data.frame(data))
    stop("data must be a data.frame")
  
  subcohort.data = data[vi==1,]
  
  #ci.method is either "standard" or "logit.transformed" 
  if(!is.element(substr(ci.method, 1,4), c("stan", "logi"))) stop("ci.method must be either 'standard' or 'logit.transformed'")

  #SEmethod is either "normal" or "boostrap"
  if(!is.element(substr(SEmethod, 1,4), c("norm", "boot"))) stop("SEmethod must be either 'normal' or 'bootstrap'")
  
  if(is.element("all", measures)) measures <- c("AUC","TPR", "FPR", "PPV", "NPV")
  
  #make sure we have a cutoff if the measures call for it
  if(any(c("TPR", "FPR", "PPV", "NPV") %in% measures) & is.na(cutoff)) stop("'cutoff' must be set in order to calculate 'FPR', 'TPR', 'PPV, 'NPV'")


    if(any(weights[vi==1] <=0)) stop("weights must be > 0.")
    if(is.element(substr(SEmethod, 1,4), c("boot"))) stop("bootstrap SE's cannot be calculated when sample weights are provided, please set SEmethod='normal'")
    subcohort = TRUE
  
  
  
  
  
  #end of checks
  cohort.size = dim(data)[1]
  
  ## get estimates via getEstimates, also calculate the bootstrap se if necessary
  
  #build data frame for getEstimates

  
  if(is.element(estimation.method, c("S", "SP", "Semi-Parametric", "semiparametric"))){

    mydata <- prepareDataSP(time, event, marker, weights, vi)  
    estRawOutput <- getEstimatesSP( data = mydata, 
                                      cutpoint = cutoff,  
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = TRUE,  
                                      cutoffN = dim(subcohort.data)[1])  


    
    
    
    
  }else if(is.element(estimation.method, c("N", "NP", "Non-Parametric", "nonparametric"))){
   # warning("Standard error calculations are not available for non-parametric estimates yet.")
    mydata <- prepareDataNP(time, event, marker, weights, vi)  
    estRawOutput <- getEstimatesNP( subcohort.data = mydata$subdata,
                                    cohort.data = mydata$cohortdata,
                                      cutpoint = cutoff,  
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = FALSE,  
                                      subcohort = TRUE)

  }else{
    
    stop("estimation.method not set correctly: it must be one of `SP` (or 'semiparametric') or 'NP' (or 'nonparametric')")
  }

  myests <- processRawOutput(estRawOutput, CImethod = ci.method, alpha)
  myests$se.cohort = estRawOutput$se.coh
  myests$fit = NULL
  myests$model.fit <- estRawOutput$fit; 
  myests$marker.cutpoint = marker.cutpoint; 
  myests$ci.method = ci.method; 
  #myests$SEmethod = SEmethod;
  myests$predict.time = predict.time; 
  myests$alpha = alpha; 
  myests$study.design = "Case-Cohort"
  myests$estimation.method <- ifelse(estimation.method=="NP", "Non-parametric", 
                                     "Semi-parametric")
  

  
  ## return the results in a nice fashion
  class(myests) <-  "SurvMTP_cch"
  myests
}
