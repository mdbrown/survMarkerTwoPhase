prepareDataSP <- function(time, event, marker, weights, vi){
  
  #probably put some checks here
  
  outData <- data.frame(time, event, marker, weights) # add vector of 1s for sample weights used later

  names(outData) <- c("xi", "di", "Y", "wi")
  outData <- outData[vi==1,] #only want the subcohort data. 
  completeCases <- complete.cases(outData)
  if(sum(completeCases) < nrow(outData)) warning(paste("NA's present, only complete cases will be used. New subcohort sample size is:", sum(completeCases)) )
  outData <- outData[completeCases, ]
  
  outData
  
}


prepareDataNP <- function(time, event, marker, weights, vi){
  
  #remove NA's from subcohort data
  tmpData <- data.frame(time, event, marker, weights) 
  tmpData <- tmpData[vi==1,]
  
  completeCases <- complete.cases(tmpData)
  if(sum(completeCases) < nrow(tmpData)) warning(paste("NA's present in subcohort data, only complete cases will be used. New subcohort sample size is:", sum(completeCases)) )
  subcohort.data <- tmpData[completeCases, ]
  
  #remove NA's from cohort data
  tmpData <- data.frame(time, event) 
  
  completeCases <- complete.cases(tmpData)
  if(sum(completeCases) < nrow(tmpData)) warning(paste("NA's present in cohort data 'time' or 'event' variables. New cohort sample size is:", sum(completeCases)) )
  cohort.data <- tmpData[completeCases, ]
  
  
  N = dim(cohort.data)[1]
  subN = dim(subcohort.data)[1]
  
  #calculate censoring weights
  psi = rep(0, subN)
  subdata.event <- subcohort.data[,2] #events in subcohort
  for( eventID in unique(subdata.event)){
     psi[is.element(subdata.event, eventID)] <- sum(is.element(cohort.data[,2], eventID))/N
  }
  
  subcohort.outData <- data.frame(cbind( subcohort.data[,c(1:3)], #time, event, marker
                                          0, subdata.event+1, psi, subcohort.data[,4]))  
  names(subcohort.outData) = c("xi","di","yi","zi","si","psi", "wi")
  
  cohort.outData        <- cohort.data[,c(1,2)]
  names(cohort.outData) <- c("xi", "di")
  out <- NULL
  out$subdata <- subcohort.outData
  out$cohortdata <- cohort.outData
  out
}

processRawOutput <- function(myests, CImethod, alpha){

if(is.null(myests$se)){ 
  myests$ci.bounds <- rbind(rep(NA, 6), rep(NA, 6));
  myests$se <- rep(NA, 6)
}else{
#calculate confidence intervals
if(substr(CImethod, 1, 4)=="stan"){
  
  myests$ci.bounds = data.frame(rbind(upper = myests$estimates - qnorm(alpha/2)*myests$se, 
                                     lower = myests$estimates - qnorm(1-alpha/2)*myests$se))
  names(myests$ci.bounds) = names(myests$estimates)
  
}else{
  #logit transform everything but the coef

  mynames = names(myests$estimates)
  myests$ci.bounds = data.frame(rbind(upper = expit(logit(myests$estimates[,mynames!="coef"]) - 
                                                     qnorm(alpha/2)*myests$se[,mynames!="coef"]/(myests$estimates[,mynames!="coef"]*(1-myests$estimates[,mynames!="coef"]))), 
                                     lower = expit(logit(myests$estimates[,mynames!="coef"]) - 
                                                     qnorm(1-alpha/2)*myests$se[,mynames!="coef"]/(myests$estimates[,mynames!="coef"]*(1-myests$estimates[,mynames!="coef"])))))
  
  if(is.element("coef", mynames)){
  myests$ci.bounds = cbind(data.frame(rbind(upper = myests$estimates[,mynames=="coef"] - qnorm(alpha/2)*myests$se[,mynames=="coef"], 
                                           lower = myests$estimates[,mynames=="coef"] - qnorm(1-alpha/2)*myests$se[,mynames=="coef"])), myests$ci.bounds)
  }
  names(myests$ci.bounds) = names(myests$estimates)
  
}
}
myests
}