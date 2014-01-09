#function to estimate summary measures for risk marker evaluation under ncc design

survMTP_ncc <- function(time, 
                        event, 
                        marker, 
                        subcoh, 
                        data, 
                        estimation.method, 
                        predict.time,
                        marker.cutpoint = median(marker)){
  
  #checks for input
  stopifnot(is.data.frame(data))

  time <- eval(substitute(time), data)
  event <- eval(substitute(event), data)
  marker <- eval(substitute(marker), data)
  s
  
  stopifnot(is.element(estimation.method, c("NP", "SP")))
  stopifnot(is.numeric(predict.time))
  stopifnot(is.numeric(marker.cutpoint))

  #ready to go
  
  #calculate weights 
  
  
  
  
  
}