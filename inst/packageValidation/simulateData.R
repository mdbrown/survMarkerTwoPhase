

#function to run a single simulation 

run_one_sim <- function(sample.design, 
                        estimation.method,
                        N = 1000,
                        beta = log(3), 
                        lam0 = 0.1, 
                        cens.perc = 0.8, 
                        ncch = 150, #for cch
                        m.match = 2, #for ncc
                        predict.time = 2, 
                        marker.cutpoint = 0){
  
  
  ##cch 
  if(sample.design=="cch"){
    
    #simulate data
    sData <- simulateData.cohort(nn = N, beta = beta, lam0 = lam0, cens.perc = cens.perc)
    sData <- simulateData.cch(sData, type = 1, ncch= 250) # draw the subsample
    
    #sample weights
    phat <- P0HAT.CCH.FUN(sData, type = 1, ncch=250)
    sData$weights <- 1/phat
    #get estimates
    ests <- survMTP.cch( time = xi, 
                         event = di, 
                         marker = Y, 
                         weights=weights, 
                         subcoh = vi, 
                         data = sData, 
                         estimation.method = estimation.method, 
                         predict.time = predict.time,
                         marker.cutpoint = marker.cutpoint)
    ests <- ests$estimates
    
  }else if(sample.design == "ncc"){
    
    #simulate data, this function returns a list with three elements
    sData <- simulateData.ncc(nn = N, beta = beta, lam0 = lam0, 
                              cens.perc = cens.perc, 
                              m.match = m.match, 
                              time.max = NULL)
    
    #build risk sets
    risk.sets <- sData[[2]]
    risk.sets <- risk.sets[,-c(1:2)] #throw away the first two columns
    
    #cohort data with subsample indicated by vi
    sData <- data.frame(sData[[1]])
    sData$id <- 1:N
    
    ests <- survMTP.ncc(time = xi, event=di, marker = yi, 
                        subcoh = vi, id = id, 
                        data = sData, 
                        risk.sets = risk.sets,
                        estimation.method = estimation.method, 
                        predict.time = predict.time, 
                        marker.cutpoint = marker.cutpoint)
    
    ests$model.fit <- NULL #dont want to save this info..too much!
    
  }else if(sample.design == "cohort"){
    stop("dont do this yet")
  }else{
    stop("sample.design not set correctly")
  }
  
  print(proc.time())
  return(ests)
  
}








######### marker simulation functions

simulateData.ncc <- function(nn=5000, 
                                  beta = beta.Ha, 
                                  lam0 = a, 
                                  cens.perc = cens.perc, 
                                  time.max = time.max, 
                                  m.match)
{
  ## ===================================================== ##
  ## Z is the matching vector; Y is the marker of interest ##
  ## ===================================================== ##
  yi = rnorm(nn)
  Zi = NULL   
  
  mu.i <- yi*beta
  
  #true survival time
  r.ti <- log(-log(runif(nn)))
  ti <-  -mu.i + r.ti
  ti <- exp(ti)/lam0
  
  #censoring time
  if(is.null(time.max)){ 
    
    if(cens.perc > 0){
      ci <- runif(nn, min = 0, max = 1)
      time.max <- uniroot(function(time.max, ti, ci, cens.perc){ mean(ti > ci*time.max) - cens.perc }, interval = c(0, 10*max(ti)), ti, ci, cens.perc)$root
      ci <- ci*time.max
    }else{
      ci <- ti+.01 # if there is no censoring, we set all the censoring times to be larger than the survival times. 
    }
  }else{
    ci = rep(time.max, nn)
  }
  
  
  
  xi = pmin(ti,ci); 
  di = 1*(ti<=ci); 
  bi = vi = di; 
  
  ind.case = (1:nn)[di==1]; 
  tivec = nivec = rep(0,sum(di)); 
  
  IND.ik = matrix(NA,nrow=sum(di),ncol=m.match+1) ## id of cases and the corresponding controls
  Iik0 = matrix(0,nrow=sum(di),ncol=nn-sum(di))  ## indicates whether xk is in the matched risk set for ti
  for(l in 1:sum(di)){
    tmpind = ind.case[l]; 
    risksetind = xi>xi[tmpind];  
    IND.ik[l,1]=tmpind;  
    ## =========================================================================== ##
    ## if matching, additional constraint of |Zi - Zl| <= a0 needs to be satisfied ##
    ## =========================================================================== ##
    
    Iik0[l,] = 1*risksetind[di==0]
    risksetind = (1:nn)[risksetind]; 
    nl = length(risksetind);
    nivec[l] = nl;
    tivec[l] = ti[tmpind];  
    ## =========================================================================== ##
    ## if riskset is empty, no control would be selected                           ##
    ## if riskset size < m.match, only select # of available for Finite Population ##
    ## =========================================================================== ##
    if(length(risksetind)>0){
      controlind = as.numeric(sample(as.character(risksetind),min(m.match,nl))); ##Finite Population
      IND.ik[l,-1] = c(controlind,rep(NA,m.match-length(controlind)))
      vi[controlind]=1
    }
    
  }
  Iik0 = Iik0[,vi[di==0]==1];  
  list("data"=cbind(xi,di,vi,yi,Zi), "Vi.k0.IND"=cbind(tivec,nivec,IND.ik),"Iik0"=Iik0) 
} 




simulateData.cohort <- 
  function(nn,
           mu = 0, 
           Sigma = 1, 
           beta = log(3), 
           lam0 = .1, cens.perc, time.max = NULL)
  {
    
    Y <- rnorm(nn, mu, Sigma)
    mu.i <- Y*beta
    
    #true survival time
    r.ti <- log(-log(runif(nn)))
    ti <-  -mu.i + r.ti
    ti <- exp(ti)/lam0
    
    #censoring time can either be set by "time.max" or by specifying cens.perc = "censoring percentage", whereas I solve for the time.max that gives a certain percentage censoring. 
    
    
    if(is.null(time.max)){ 
      
      if(cens.perc > 0){
        ci <- runif(nn, min = 0, max = 1)
        time.max <- uniroot(function(time.max, ti, ci, cens.perc){ mean(ti > ci*time.max) - cens.perc }, interval = c(0, 10*max(ti)), ti, ci, cens.perc)$root
        
        ci <- ci*time.max
      }else{
        ci <- ti+.01 # if there is no censoring, we set all the censoring times to be larger than the survival times. 
      }
      
    }else{
      
      #ci <- runif(nn, min = 0, max = time.max)
      ci = rep(time.max, nn)
    }
    
    #observed marker is the min of ti and ci        
    xi <- pmin(ti, ci)
    # failure indicator
    di <- ifelse( ti == xi, 1, 0)
    
    #xi is the min of ti and ci
    #di is the indicator for failure, 1 for failure, 0 for censored
    #Y is the marker values
    
    result <- as.data.frame(cbind(xi, di, Y)) 
    names(result) = c( "xi", "di", "Y")
    return(result)
  }



###
### cch
###

## type = 1: sample all case, and a random sample from the full cohort 
## type = 2: sample without replacement among di=0 and di=1  

simulateData.cch <- function(data0, ncch  = NULL,
                        ncch0 = NULL,
                        ncch1 = NULL,
                        type = 2)
{ 
  nn <- nrow(data0)
  vi <- rep(0,nn)
  di <- data0$di
  
  if (type==1) {
    
    vi <- di # choose all cases
    
    junk.ind <- sample(1:nn,ncch) #sample controls
    
    vi[junk.ind] <- 1; 
    
  } else {
    
    case.ind    <- (1:nn)[di==1]
    control.ind <- (1:nn)[di==0]
    
    ncase    <- length(case.ind)
    ncontrol <- length(control.ind)
    
    temp.ncch0 <- min(ncontrol, ncch0) #min of what we have and how many we want
    temp.ncch1 <- min(ncase, ncch1)
    
    junk.ind0 <- sample(ncontrol, temp.ncch0, replace = FALSE) #sample controls
    vi[control.ind[junk.ind0]] <- 1 
    
    junk.ind1 <- sample(ncase, temp.ncch1, replace = FALSE)    #sample cases
    vi[case.ind[junk.ind1]] <- 1
  }
  
  cohortdata    <- data0
  cohortdata$vi <- vi
  
  cohortdata  
}


P0HAT.CCH.FUN<-function(cohortdata,type,ncch, ncch0,ncch1) {
  
  nn <- dim(cohortdata)[1]
  di <- cohortdata$di
  
  if(type==1){ 
    p0hat <- di+ (1- di)*ncch/nn
  } else{
    p0hat <- di*ncch1/sum(di) + (1-di)*ncch0/(nn-sum(di)) 
  }
  
  p0hat
}  
