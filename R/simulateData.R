
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