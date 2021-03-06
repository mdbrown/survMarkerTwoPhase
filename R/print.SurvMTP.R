#' Print an object of class SurvAM--output from survEstMeasures
#' 
#'  
#' @param x,output from the function survEstMeasures
#' @return NULL


print.SurvMTP_cch <- function(x, ...){
  #methods to print an object of class "SurvAM"
  # want to display the estimates, standard errors, and confidence intervals
  # just like coxph
  
  #x is a list with elements 'estimates', 'se', 'CIbounds', 'cutoff', 'CImethod', 'SEmethod', and 'predict.time'
  cat(paste("\n",x$estimation.method, "estimates under", x$study.design, "study design.\n"))
  cat("\n")

  mynames = names(x$estimates)
  
  if(any(mynames %in% c("FPR", "TPR" , "NPV" , "PPV"))){
    whitespace <- rep(" ",4)
  }else{
    whitespace <- rep(" ", 4)
  }
  
  cat(whitespace, #"estimate")
            paste("estimate     se      lower ",
            1-x$alpha, "  upper ",
            1-x$alpha, "\n",sep = ""))
  
  cat("\n")

  for(i in 1:length(mynames)){
    if(mynames[i] %in% c("FPR", "TPR" , "NPV" , "PPV") ) mynames[i] = paste(mynames[i], "(c)", sep = "")
    cat(paste(sprintf("%-6s", mynames[i]), 
              sprintf("%10.3f", round(x$estimate[i], 3)), 
              sprintf("%10.3f ", round(x$se[i], 3)), 
              sprintf("%13.3f ", round(x$ci.bounds[2,i], 3)), 
              sprintf("%11.3f ", round(x$ci.bounds[1,i], 3)),
               sep = "")); 
    cat("\n")
    
  }
  
  cat("\n")
  
if(any(mynames %in% c("FPR(c)", "TPR(c)" , "NPV(c)" , "PPV(c)"))) cat(" marker cutpoint: c =", x$marker.cutpoint, "\n")
cat("\n")
  
  
}




print.SurvMTP_ncc <- function(x, ...){
  #methods to print an object of class "SurvAM"
  # want to display the estimates, standard errors, and confidence intervals
  # just like coxph
  
  #x is a list with elements 'estimates', 'se', 'CIbounds', 'cutoff', 'CImethod', 'SEmethod', and 'predict.time'
  cat(paste("\n",x$estimation.method, "estimates under", x$study.design, "study design.\n"))
  cat("\n")
  
  mynames = names(x$estimates)
  
  if(any(mynames %in% c("FPR", "TPR" , "NPV" , "PPV"))){
    whitespace <- rep(" ",4)
  }else{
    whitespace <- rep(" ", 4)
  }
  
  cat(whitespace, #"estimate",
  paste("estimate     se      lower ",
  1-x$alpha, "  upper ",
  1-x$alpha, "\n",sep = ""))
  
  cat("\n")
  
  for(i in 1:length(mynames)){
    if(mynames[i] %in% c("FPR", "TPR" , "NPV" , "PPV") ) mynames[i] = paste(mynames[i], "(c)", sep = "")
    cat(paste(sprintf("%-6s", mynames[i]), 
              sprintf("%10.3f", round(x$estimate[i], 3)), 
              sprintf("%10.3f ", round(x$se[i], 3)), 
              sprintf("%13.3f ", round(x$ci.bounds[2,i], 3)), 
              sprintf("%11.3f ", round(x$ci.bounds[1,i], 3)),
              sep = "")); 
    cat("\n")
    
  }
  
  cat("\n")
  
  if(any(mynames %in% c("FPR(c)", "TPR(c)" , "NPV(c)" , "PPV(c)"))) cat(" marker cutpoint: c =", round(x$marker.cutpoint,3), "\n")
  cat("\n")
  
  
}
