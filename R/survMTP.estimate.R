 
survMTP.estimate <- function(time, event, marker, weights, 
                             subcohort.data,
                             cohort.data = NULL, 
                             cohort.size, 
                             estimation.method = "NP", 
                             study.design = "Case-Cohort", 
                             predict.time, 
                             measures = c('all'), 
                             cutpoint = median(marker)
                           ){
  cutoff <- cutpoint
  cutoff.type = "none"; #cutoffN = 100 #functionality to be added later
  SEmethod ="normal"; #bootstraps = 10
  CImethod = "logit.transformed"
  alpha=0.05
  
  #put checks here
  if(length(cutoff)==0) cutoff = NA;  
  
  #data.frame checks
  if(!is.data.frame(subcohort.data))
    stop("subcohort.data must be a data.frame")
  
  if(estimation.method=="NP" & is.null(cohort.data))
    stop("cohort.data must be supplied when using estimation method NP")
  
  if(!is.null(cohort.data)){
    if(!is.data.frame(cohort.data))
       stop("cohort.data must be a data.frame")
  }
  
  
  #CImethod is either "standard" or "logit.transformed" 
  if(!is.element(substr(CImethod, 1,4), c("stan", "logi"))) stop("CImethod must be either 'standard' or 'logit.transformed'")

  #SEmethod is either "normal" or "boostrap"
  if(!is.element(substr(SEmethod, 1,4), c("norm", "boot"))) stop("SEmethod must be either 'normal' or 'bootstrap'")
  
  if(is.element("all", measures)) measures <- c("AUC","TPR", "FPR", "PPV", "NPV")
  
  #make sure we have a cutoff if the measures call for it
  if(any(c("TPR", "FPR", "PPV", "NPV") %in% measures) & is.na(cutoff)) stop("'cutoff' must be set in order to calculate 'FPR', 'TPR', 'PPV, 'NPV'")

  if(is.null(weights)){ 
    stop("Must specify weights!")
   # weights = rep(1, N)
   # subcohort = FALSE
  }else{
    if(any(subcohort.data[,weights] <=0)) stop("weights must be > 0.")
    if(is.element(substr(SEmethod, 1,4), c("boot"))) stop("bootstrap SE's cannot be calculated when sample weights are provided, please set SEmethod='normal'")
    subcohort = TRUE
  }
  
  #end of checks
  
  
  ## get estimates via getEstimates, also calculate the bootstrap se if necessary
  
  #build data frame for getEstimates

  
  if(is.element(estimation.method, c("S", "SP", "Semi-Parametric", "semiparametric"))){

    mydata <- prepareDataSP(time, event, marker, weights, 
                            subcohort.data)  
    estRawOutput <- getEstimatesSP( data = mydata, 
                                      cutpoint = cutpoint,  
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = FALSE,  
                                      cutoffN = dim(subcohort.data)[1])  
    myests <- NULL
    myests$estimates <- estRawOutput$estimates
    
    names(myests$estimates) = c("coef", measures)
    #names(myests$se) = c("coef", measures)
    myests$model.fit <- estRawOutput$fit; 
    myests$cutpoint = cutoff; 
    #myests$CImethod = CImethod; 
    #myests$SEmethod = SEmethod;
    myests$predict.time = predict.time; 
    #myests$alpha = alpha; 
    myests$study.design = study.design
    myests$estimation.method = "Semi-parametric"
  }else if(is.element(estimation.method, c("N", "NP", "Non-Parametric", "nonparametric"))){

    mydata <- prepareDataNP(time, event, marker, weights, 
                            subcohort.data, 
                            cohort.data)  
    estRawOutput <- getEstimatesNP( subcohort.data = mydata$subdata,
                                    cohort.data = mydata$cohortdata,
                                      cutpoint = cutpoint,  
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = FALSE,  
                                      subcohort = TRUE)
    myests <- estRawOutput
    
    names(myests$estimates) = c( measures)
    #names(myests$se) = c("coef", measures)
    #myests$model.fit <- myests$fit; 
    myests$cutpoint = cutoff; 
    #myests$CImethod = CImethod; 
    #myests$SEmethod = SEmethod;
    myests$predict.time = predict.time; 
    #myests$alpha = alpha; 
    myests$study.design = study.design
    myests$estimation.method = "Non-parametric"
  }else{
    
    stop("estimation.method not set correctly: it must be one of `SP` (or 'semiparametric') or 'NP' (or 'nonparametric')")
  }
  
  
  

  
  ## return the results in a nice fashion
  class(myests) <-  "SurvMTP"
  myests
}
