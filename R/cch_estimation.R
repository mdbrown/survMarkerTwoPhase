#' get Risk estimates and standard errors, called from SurvRM

#'@param data data.frame of xi (surv time), di (event ind), and Y (marker)
#'@param measures character vector of measures wanted, can be a subset of c("AUC", "beta", "TPR", "FPR", "PPV", "NPV")
#'@param predict.time numeric prediction time to evaluate measures at
#'@param CalVar Logical, should standard errors be calculated?
#'@return a list consisting of 'est' with point estimates in order of 'measures' and 
#'        (if CalVar == TRUE), 'sd' with standard errors for estimates. 


getEstimatesSP <- function(data, 
                           cutpoint,  
                           measures,
                           predict.time,
                           CalVar=TRUE, cutoff.type = "none", cutoffN = 100, subcohort=FALSE)
{  
  
  #  browser()
  
  N = nrow(data)
  data$vi = 1; #data$wi = 1
  cutoff <- cutpoint
  
  #junk = GetRTdata(data, predict.time)  ## data.RT and data are sorted by Y 
  ###
  
  ####
  ## Est.Sy
  
  
  data.s <- data[order(-data$xi),] #data sorted by stime
  
  fit<-coxph( Surv(data$xi,data$di) ~ data$Y, weights = data$wi, method = "breslow")   #original data, sorting doesn't matter. 
  betahat<-fit$coef                                               
  
  
  r.riskset <-data.s$wi/cumsum(exp(data.s$Y*betahat )*data.s$wi)
  Lambda0t  <- sum(r.riskset[(data.s$xi <= predict.time)&(data.s$di == 1)])
  
  
  
  linearY <-  data$Y*betahat            #linearY is sorted by original data
  Sy <- exp(-Lambda0t*exp(linearY))
  
  
  ####
  
  
  ooo = order( data$Y)
  data.RT <- cbind(data$Y, Sy, data$wi)[ooo,] 
  
  Fck <- sum.I(data.RT[,1], ">=", data.RT[,1], data.RT[,3])/sum(data.RT[,3])
  
  data.RT <- cbind(data.RT[,-c(3)],Fck) 
  
  ###

  cutoff.type = "yes"
  if(cutoff.type != "none"){
    
    #if we use cutoffs, not used in the released version of package
    cutoffN <- min(N, 100)
    cutoffs <- unique(sort(c( cutpoint, quantile(linearY, (1:cutoffN/cutoffN), type =1, na.rm = TRUE))))
    
    cutpos = sum.I(cutoffs,">=", linearY[ooo])
    subdata.RT = data.RT[cutpos, ]
    
    RT.out = EstRTall(subdata.RT) 
    
    AUC    <- RT.out[[2]]
    RT.out <- RT.out[[1]]
    
    
  }else{
    
    RT.out = EstRTall(data.RT) 
    AUC       <- RT.out[[2]]  
    RT.out    <- RT.out[[1]]
    
  }
  
  
  if (length(measures[measures!="AUC"])>0) {
    typey = measures[measures!="AUC"]; 
    typex = rep("cutoff", length(typey))
    vp = rep(cutoff, length(typey));
    
    tmpind <- with(RT.out, which.min(abs(cutoff-cutpoint)))
    all.measures = RT.out[tmpind,]
    
    
    RTvp.out  = all.measures[, measures[measures!="AUC"]]
    
  }else{
    RTvp.out = NULL; typex = typey = vp = NULL; 
  }
  
  
  if("AUC" %in% measures) est = unlist(c(betahat, AUC, RTvp.out ))  else est = unlist(c(betahat, RTvp.out))
  
  est <- data.frame(t(est)); 
  names(est) = c("coef", measures)
  
  if (CalVar)  {
    
    subdata = cbind(data[ooo,],data.RT[,c(2)], linearY[ooo])
    names(subdata)=c("times","status","y","wi","vi","Sy","linearY")
    if(cutoff.type=="none") cutoffs = NA
    jjunk = Est.Wexp.cpp(subdata, 
                         N,
                         RT.out,
                         predict.time,
                         vp,
                         typex,
                         typey, 
                         resid(fit, "score")[ooo], 
                         fit$var, 
                         cutoffs = cutoffs)
    
    Wexp = data.frame(cbind(jjunk$Wexp.beta,jjunk$Wexp.AUC,jjunk$Wexp.vp))

    tmpse = Est.Var.CCH.trueweights(N,Wexp,subdata,subdata$status, subcohort=TRUE) 
    
    se <- sqrt(tmpse$cch.variance)
    se.coh <- sqrt(tmpse$cohort.variance)
    
    se <- data.frame(t(se))
    se.coh<- data.frame(t(se.coh))
    names(se) = c("coef", measures)
    names(se.coh) = c("coef", measures)
    list(estimates = est, se = se, fit = fit, se.coh = se.coh) 
    
    
  } else {list(estimates = data.frame(est), fit = fit)}
}







getEstimatesNP <- function(subcohort.data, 
                           cohort.data,
                           cutpoint,  
                           measures,
                           predict.time,
                           CalVar = TRUE, subcohort=FALSE)
{
  
  
  
  
  
  N = nrow(cohort.data); ## cohort size
  data0 <- cohort.data #unsorted data
  data  <- subcohort.data
  data = data [order(data$yi),] ## sorted by marker
  
  ck = data$yi;     
  wgtk = data$wi;  
  xk = data$xi; 
  sk = data$si; 
  psk=data$psi; 
  
  c0 <- cutpoint; 
  t0 <- predict.time
  
  
  nc = length(ck); # sampled size 
  ind.ck = (1:nc)[order(ck)]
  CWk = WGT.FUN(data[,c(1,2)],data0, t0=t0)  # use full cohort xi di to calculate censoring weights
  
  St0.Fck = sum.I(ck,">=",ck,wgtk*CWk*(xk >= t0))/sum(CWk*wgtk)
  Ft0.Fck = sum.I(ck,">=",ck,wgtk*CWk*(xk <  t0))/sum(CWk*wgtk)
  Fck = sum.I(ck,">=",ck,wgtk*CWk)/sum(CWk*wgtk)
  St0 = max(St0.Fck)            ## St0     = P(T> t0)
  Ft0 = 1-St0                   ## Ft0     = P(T<=t0)
  FPR.ck= (St0-St0.Fck)/St0     ## P(Y> ck|T> t0)
  TPR.ck= (Ft0-Ft0.Fck)/Ft0     ## P(Y> ck|T<=t0)
  NPV.ck= St0.Fck/Fck           ## P(T> t0|Y<=ck)
  PPV.ck= (Ft0-Ft0.Fck)/(1-Fck) ## P(T<=t0|Y> ck)
  AUC = sum(TPR.ck*(FPR.ck-c(FPR.ck[-1],0)))
  
  ## acc.ck: accuracy estimates at c
  nm.acc = c("FPR","TPR","NPV","PPV"); 
  acc.ck = data.frame("cutoff"=ck,"FPR"=FPR.ck,"TPR"=TPR.ck,"NPV"=NPV.ck, "PPV"=PPV.ck)    
  
  
  if(!is.null(c0)){
    tmpind.c0 = sum.I(c0, ">=", ck); acc.c0 = acc.ck[tmpind.c0,]; F.c0 = Fck[tmpind.c0]
  }else{ acc.c0 = acc.ck }
  est= data.frame("AUC" = AUC, acc.c0[-1])
  
  est  = est[,measures]
  #est= list("AUC" = AUC, "ACC.u0"=acc.uk[tmpind.u0,-c(1,ind0)],"ACC.c0"=acc.c0, "ACC.all" = acc.ck) ##output all cutoff
  if(!CalVar){
    list("estimates" = est)
  }else{
    ###### Variance calculation below ########
    Phi=Phi.C.new.FUN(xk=data$xi,dk=data$di, Ti=data0$xi, Di=data0$di, t0 = t0)
    
    ## doing u0 and c0 together
    c.u0 = NULL; acc.u0.temp=NULL; F.c0.b =F.c0
    
    CC = c(c.u0,c0); #nu0 = length(u0); 
    acc.c0.b = rbind(acc.u0.temp,acc.c0); 
    
    U.ACC.c0.tmp = as.list(1:4); 
    U.ACC.c0 =Wexp.c0= as.list(1:4); 
    names(U.ACC.c0.tmp) = names(U.ACC.c0)=nm.acc; n.acc.c=length(U.ACC.c0)
    
    
    I.ck.c0 = 1*(ck>=VTM(CC,nc)); 
    U.ACC.c0.tmp$FPR =  (xk >  t0)*(I.ck.c0-  VTM(acc.c0.b$FPR,nc))/St0      ## exp for FPRhat(c)-FPR(c)
    U.ACC.c0.tmp$TPR =  (xk <= t0)*(I.ck.c0-  VTM(acc.c0.b$TPR,nc))/(1-St0)  ## exp for TPRhat(c)-TPR(c)
    U.ACC.c0.tmp$NPV =  (1-I.ck.c0)*(1*(xk> t0)-VTM(acc.c0.b$NPV,nc))/VTM(F.c0.b,nc)
    U.ACC.c0.tmp$PPV =    I.ck.c0*(1*(xk<=t0)-VTM(acc.c0.b$PPV,nc))/(1-VTM(F.c0.b,nc))
    
    U.AUC = (xk<=t0)/(1-St0)*(1-FPR.ck-AUC)+(xk>t0)/St0*(TPR.ck-AUC)
    
    Wexp.np.AUC = CWk*U.AUC+Phi%*%(wgtk*CWk*U.AUC)/sum(wgtk)
    se.auc = Est.Var.CCH.trueweights(N,Wexp.np.AUC,data,data$si, subcohort=TRUE)
    se.u0 = NULL
    
    se.c0 = NULL
    se.coh = NULL
    
    if (!is.null(c0)) {
      npc = length(U.ACC.c0)  
      for(kk in 1:npc){ 
        
        U.ACC.c0[[kk]] = U.ACC.c0.tmp[[kk]]
        Wexp.c0[[kk]] = CWk*U.ACC.c0[[kk]]+Phi%*%(wgtk*CWk*U.ACC.c0[[kk]])/sum(wgtk)
        tmp.se <- Est.Var.CCH.trueweights(N,data.frame(Wexp.c0[[kk]]),data,data$si, subcohort=subcohort)
        se.c0  = c(se.c0,tmp.se$cch.variance)
        se.coh = c(se.coh, tmp.se$cohort.variance) 
      }
      se.c0 = data.frame(sqrt(matrix(se.c0,nrow=length(c0))))
      se.coh = data.frame(sqrt(matrix(se.coh, nrow = length(c0))))
      names(se.c0) = names(se.coh) = nm.acc    
    }
    
    se <- data.frame("AUC" = sqrt(se.auc$cch.variance), se.c0)
    se.coh <- data.frame("AUC" = sqrt(se.auc$cohort.variance), se.coh)
    
    se <- se[,measures]
    se.coh <- se.coh[,measures]
    list("estimates" = est,"se" =se, "se.coh" = se.coh) 
  }	
  
  
}
