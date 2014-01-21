#function to estimate summary measures for risk marker evaluation under ncc design

survMTP.ncc <- function(time, 
                        event, 
                        marker, 
                        subcoh, 
                        id, 
                        data,
                        risk.sets, 
                        estimation.method = "NP", 
                        predict.time,
                        alpha = 0.05, 
                        Npert = 500, 
                        marker.cutpoint = 'median' )
  {
  

  #checks for input
  stopifnot(is.data.frame(data))
  
  time <- eval(substitute(time), data)
  event <- 1*eval(substitute(event), data)
  marker <- eval(substitute(marker), data)
  vi <- 1*eval(substitute(subcoh), data)
  id <- eval(substitute(id), data)
  
  stopifnot(is.element(estimation.method, c("NP", "SP")))
  stopifnot(is.numeric(predict.time))
  
  if(marker.cutpoint=='median') marker.cutpoint =  median(eval(substitute(marker), data))
  stopifnot(is.numeric(marker.cutpoint))

  # remove missing data if time, event, id, subcohort is missing

  completeIND <- complete.cases(time, event, vi, id)
  
  if(!all(completeIND)){ 

    stop(paste('There are missing values in the 
                   time, event, subcoh or id variables.
                   Please remove these from "data" and re-run' ))
  
#  time  <- time[completeIND]
#  event <- event[completeIND]
#  vi <- vi[completeIND]
#  id <- id[completeIND]
#  marker <- marker[completeIND]
  }
  
  #remove observations with missing marker value if it is in the subcohort
  
  completeIND_sub <- complete.cases(marker[vi==1])
  
  if(!all(completeIND_sub)){ 
    
    stop(paste('There are missing marker values in the subcohort. Please remove and re-run'))
    #ids to remove
   # rmID <- id[vi==1][!completeIND_sub]
    
   # time  <- time[!is.element(id, rmID)]
  #  event <- event[!is.element(id, rmID)]
  #  vi <- vi[!is.element(id, rmID)]
  #  marker <- marker[!is.element(id, rmID)]
  #  id <- id[!is.element(id, rmID)]
  }
  
  #\checks
  
  #build the data structures needed to call the necessary functions
  
  #calculate risk set sizes
  case.IDs <- id[event==1]
  case.times <- time[event==1]

  risk.set.N <- sapply( case.times, function(x) sum(time > x))
  
  #so i can match the risk.set.N and time up with the already provided risk.sets
  ind <- match(case.IDs, risk.sets[,1])
  risk.sets <- data.frame(case.times[ind], risk.set.N[ind], risk.sets)
  
  ## ixk indicator matrix, one row for each case, one column for each control  
  ## indicates whether control k is in the matched risk set for case i 
  N <- nrow(data)
  Iik <- matrix(0,nrow=sum(event),ncol= sum(event==0 & vi==1) )  
  control.times <-time[event==0 & vi ==1]
  
  for( i in 1:nrow(Iik)){
    
    Iik[i, ] = 1*c(case.times[i] <  control.times)
  }
  
  inputData <- data.frame(xi = time, di = event, vi = vi, yi = marker)
  
  
  if(estimation.method=="NP"){
    

    result <- NCC_EstandSE_NP(data = list(inputData, risk.sets, Iik), 
                             cutpoint = marker.cutpoint,
                             predict.time = predict.time, 
                             nmatch = ncol(risk.sets)-3, 
                             B0 = Npert)

  }else if(is.element(estimation.method, c("S", "SP", "Semi-Parametric", "semiparametric"))){

    result <- NCC_EstandSE_SP(data = list(inputData, risk.sets, Iik), 
                             cutpoint = marker.cutpoint,
                             predict.time = predict.time, 
                              nmatch = ncol(risk.sets)-3, 
                              B0 = Npert)

    
  }else{
    stop("estimation.method not set correctly: must be either 'NP' for non-parametric or 'SP' for semi-parametric")
  }

  myests <- processRawOutput(result, CImethod = "standard", alpha)
  myests$fit = NULL
  myests$model.fit <- result$fit; 
  myests$marker.cutpoint = marker.cutpoint; 
  #myests$CImethod = CImethod; 
  #myests$SEmethod = SEmethod;
  myests$predict.time = predict.time; 
  myests$alpha = alpha; 
  myests$study.design = "Nested Case-Control"
  myests$estimation.method <- ifelse(estimation.method=="NP", "Non-parametric", 
                                                              "Semi-parametric")
  
  

  class(myests) <-  "SurvMTP_ncc"
  
  myests
  
}