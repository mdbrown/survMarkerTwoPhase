prepareDataSP <- function(time, event, marker, weights, 
                          subcohort.data){
  
  #probably put some checks here
  
  outData <- as.data.frame(subcohort.data[,c(time, event, marker, weights)]) # add vector of 1s for sample weights used later

  names(outData) <- c("xi", "di", "Y", "wi")
  completeCases <- complete.cases(outData)
  if(sum(completeCases) < nrow(outData)) warning(paste("NA's present, only complete cases will be used. New subcohort sample size is:", sum(completeCases)) )
  outData <- outData[completeCases, ]
  
  outData
  
}


prepareDataNP <- function(time, event, marker, weights, 
                          subcohort.data, 
                          cohort.data){
  
  #remove NA's from subcohort data
  tmpData <- as.data.frame(subcohort.data[,c(time, event, marker, weights)]) # add vector of 1s for sample weights used later
  
  completeCases <- complete.cases(tmpData)
  if(sum(completeCases) < nrow(tmpData)) warning(paste("NA's present in subcohort data, only complete cases will be used. New subcohort sample size is:", sum(completeCases)) )
  subcohort.data <- tmpData[completeCases, ]
  
  #remove NA's from cohort data
  tmpData <- as.data.frame(cohort.data[,c(time, event)]) 
  
  completeCases <- complete.cases(tmpData)
  if(sum(completeCases) < nrow(tmpData)) warning(paste("NA's present in cohort data, only complete cases will be used. New cohort sample size is:", sum(completeCases)) )
  cohort.data <- tmpData[completeCases, ]
  
  
  N = dim(cohort.data)[1]
  subN = dim(subcohort.data)[1]
  
  #calculate censoring weights
  psi = rep(0, subN)
  subdata.event <- subcohort.data[,event]
  for( eventID in unique(subdata.event)){
     psi[is.element(subdata.event, eventID)] <- sum(is.element(cohort.data[,event], eventID))/N
  }
  
  subcohort.outData <- data.frame(cbind( subcohort.data[,c(time, event, marker)],
                                          0, subdata.event+1, psi, subcohort.data[,weights]))  
  names(subcohort.outData) = c("xi","di","yi","zi","si","psi", "wi")
  
  cohort.outData        <- cohort.data[,c(time, event)]
  names(cohort.outData) <- c("xi", "di")
  out <- NULL
  out$subdata <- subcohort.outData
  out$cohortdata <- cohort.outData
  out
}

processRawOutput <- function(myests, CImethod, alpha){


#calculate confidence intervals
if(substr(CImethod, 1, 4)=="stan"){
  
  myests$CIbounds = data.frame(rbind(upper = myests$estimates - qnorm(alpha/2)*myests$se, 
                                     lower = myests$estimates - qnorm(1-alpha/2)*myests$se))
  names(myests$CIbounds) = names(myests$estimates)
  
}else{
  #logit transform everything but the coef

  mynames = names(myests$estimates)
  myests$CIbounds = data.frame(rbind(upper = expit(logit(myests$estimates[,mynames!="coef"]) - 
                                                     qnorm(alpha/2)*myests$se[,mynames!="coef"]/(myests$estimates[,mynames!="coef"]*(1-myests$estimates[,mynames!="coef"]))), 
                                     lower = expit(logit(myests$estimates[,mynames!="coef"]) - 
                                                     qnorm(1-alpha/2)*myests$se[,mynames!="coef"]/(myests$estimates[,mynames!="coef"]*(1-myests$estimates[,mynames!="coef"])))))
  
  if(is.element("coef", mynames)){
  myests$CIbounds = cbind(data.frame(rbind(upper = myests$estimates[,mynames=="coef"] - qnorm(alpha/2)*myests$se[,mynames=="coef"], 
                                           lower = myests$estimates[,mynames=="coef"] - qnorm(1-alpha/2)*myests$se[,mynames=="coef"])), myests$CIbounds)
  }
  names(myests$CIbounds) = names(myests$estimates)
  
}

myests$model.fit <- myests$fit; 

myests
}