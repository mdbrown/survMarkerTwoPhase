## data: xi, di, vi, yi; 
NCC.Cox.explicit.PTB.s.y <- function(data,V.IND,Iik0=Iik0.mat,wgtk.ptb=NULL,B0=500,Zmatch.ind,t0, m.match)
{
  xi = data[,1]; 
  di = data[,2]; 
  vi = data[,3]; 
  yi = data[,-(1:3),drop=F];
  yi=as.matrix(yi); 
  n.t0 = length(t0)
  
  NN = length(xi); 
  nv = sum(vi); 
  tj = V.IND[,1]; 
  pi.tj = V.IND[,2]/NN
  wgtk = 1/FNCC.WGT.FUN(data,V.IND,Iik0=Iik0.mat,Zmatch.ind, m.match = m.match);#print(sum(wgtk))
  
  xk = xi[vi==1]; 
  dk = di[vi==1]; 
  yk = yi[vi==1,,drop=F]; 
  idk = (1:NN)[vi==1];  
  
  #wgtk.ptb = PtbNCC.WGT.FUN(data,V.IND,Iik0,Zmatch.ind=Zmatch.ind,B0=B0); wgtk.ptb = wgtk.ptb[match(idk,wgtk.ptb[,1]),-1]
  wgtk.ptb = PtbNCC.WGT.FUN(data,V.IND,Iik0,B0=B0, m.match = m.match); 
  
  wgtk.ptb = wgtk.ptb[match(idk,wgtk.ptb[,1]),-1]
  
  fit = coxph(Surv(xk,dk)~yk + cluster(seq_along(yk)),weights=wgtk) #,robust=T); 
  betahat=fit$coef; ebyk = c(exp(yk%*%betahat));
  pi.tj = PI.k.FUN(tj,ebyk,xk,yk,wgtk)/sum(wgtk)
  LamCox.t0 = sum.I(t0, ">=", tj, 1/pi.tj)/sum(wgtk); gammahat = c(log(LamCox.t0),betahat)
  gammaptb = PTB.Gamma.FUN.explicit(wgtk.ptb,cbind(xk,dk,wgtk,yk),gammahat=gammahat, t0 = t0) 
  gamma.sd = apply(gammaptb,1,sd)
  #z95.suplogLam = quantile(apply(abs(gammaptb[1:n.t0,]-gammahat[1:n.t0])/gamma.sd[1:n.t0],2,max),0.95,na.rm=T)
  
  y0.cut=data[,-(1:3)][data[,3]==1]
  y0.cut=as.matrix(y0.cut); n.y0 = nrow(y0.cut); ## n.y x n.y0
  logLamptb.t0.y0 = gammaptb[rep(1:n.t0,n.y0),] + (y0.cut%*%gammaptb[-(1:n.t0),])[rep(1:n.y0,rep(n.t0,n.y0)),]
  logLamhat.t0.y0 = gammahat[rep(1:n.t0,n.y0) ] + (y0.cut%*%gammahat[-(1:n.t0)] )[rep(1:n.y0,rep(n.t0,n.y0))] 
  logLamhat.t0.y0.sd = apply(logLamptb.t0.y0,1,sd)
  yk.old <- yk
  yk = ebyk
  Fyk = sum.I(yk, ">=", yk, wgtk)/sum(wgtk) 
  Fyk.ptb = matrix(0,length(yk),B0)
  
  for (i in 1:B0) {Fyk.ptb[,i] =sum.I(yk, ">=", yk, wgtk.ptb[,i])/sum(wgtk.ptb[,i])}
  #Fyk.ptb = apply(wgtk.ptb,2,myfyk(yk,wgtk))
  
  list(ck =yk, 
       CondSck = exp(-exp(logLamhat.t0.y0)),
       CondSck.ptb = exp(-exp(logLamptb.t0.y0)), 
       beta.est=betahat,
       beta.sd =gamma.sd[-c(1:n.t0)],
       Fck = Fyk,
       Fck.ptb = Fyk.ptb, 
       yk = yk.old, 
       fit = fit)
}

FNCC.WGT.FUN <- function(data,V.IND, Iik0, Zmatch.ind=NULL, m.match)
{
  xi = data[,1];  di = data[,2]; vi = data[,3]; nv = sum(vi); n1 = sum(di)
  xk = xi[vi==1]; dk = di[vi==1]; phatk = dk; ind0 = dk == 0; nv0 = sum(ind0)
  ## tl: event time for the cases; nl: (matched) risk set size for tl
  tl = V.IND[,1]; n.Rl = V.IND[,2]; 
  ooo = order(tl); tl=tl[ooo]; n.Rl=n.Rl[ooo]
  if(is.null(Zmatch.ind))
  {
    pl = c(1,cumprod(1-m.match/(n.Rl)))
    phatk[ind0]=1-pl[sum.I(xk[ind0],">",tl)+1]
  }else{
    junk = Iik0*log(1-m.match/pmax(n.Rl,m.match)); junk[is.na(junk)] = 0
    phatk[ind0]=1-exp(apply(junk,2,sum))
  }
  phatk
}

PtbNCC.WGT.FUN <- function(data, V.IND, Iik0, Zmatch=NULL,B0=500,wgtk = NULL, m.match)
{ 

  ## V.IND: 1st col = ti; 2nd col = n(ti);  3:(nv0+2): indicate xk is the control of ti ##
  ## Iik0: indicates the matched risk set of ti 
  ti = V.IND[,1]; n.Ri = V.IND[,2]; id.ii = V.IND[,3]; V.IND = V.IND[,-(1:3)]; n1 = length(ti)
  dj.i = xj.i = matrix(NA,ncol=m.match,nrow=n1)
  
  V.IND <- as.matrix(V.IND)
  #xj.i[!is.na(c(V.IND))] = data[c(V.IND),1] 
  xj.i= data[c(V.IND),1] ;
  xj.i=xj.i[!is.na(xj.i)]
  #dj.i[!is.na(c(V.IND))] = data[c(V.IND),2] 
  dj.i= data[c(V.IND),2] ; dj.i=dj.i[!is.na(xj.i)]
  ind.control = c(dj.i!=1) ## try to identify those dj.i = 1 and sampled as control 
  id.0j.all = sort(unique(na.omit(V.IND[ind.control]))); id.all = sort(unique(c(id.ii,id.0j.all)))
  x0j = c(xj.i)[match(id.0j.all,c(V.IND))]; nsub = length(id.all)
  
  wptb.mat = matrix(0,ncol=B0,nrow=nsub); 
  for(bb in 1:B0)
  {
    Vjptb = rep(0,nsub); Bii = rexp(n1); B0mat = matrix(rexp(n1*m.match),ncol=m.match) 
    V0jptb = 1 - tapply(c(1-B0mat[ind.control]), c(V.IND[ind.control]), prod,na.rm=T)
    Vjptb[match(id.0j.all,id.all)] = V0jptb	
    Vjptb[match(id.ii,    id.all)] = Bii	
    n.Ri = ifelse(n.Ri==0,1,n.Ri)
    dLam.ti = apply(B0mat*!is.na(V.IND), 1, sum)/(n.Ri);
    Lam0j.ptb = apply(dLam.ti*Iik0,2,sum)
    pjptb = rep(1,nsub); pjptb[match(id.0j.all,id.all)] = 1-exp(-Lam0j.ptb)
    wptb.mat[,bb] = pmax(0,Vjptb/pjptb);  
  }
  if(is.null(wgtk))
  {
    return(cbind(id.all,wptb.mat)[order(id.all),])
  }else{
    return(list("wptb"=cbind(id.all,wptb.mat)[order(id.all),],"wptb.naive"=wgtk*matrix(rexp(nsub*B0),nrow=nsub)))
  }
}

PI.k.FUN <- function(tt,ebyi,xi,yi,wgt.ptb=NULL,k0=0)
{
  out = ebyi; yi=as.matrix(yi); py = ncol(yi); 
  if(!is.null(wgt.ptb))
  { 
    wgt.ptb=as.matrix(wgt.ptb); pv=ncol(wgt.ptb); 
    if(k0==1){yi=yi[,rep(1:py,pv),drop=F]; out=out*wgt.ptb[,rep(1:pv,rep(py,pv))] # v1*z1, v1*z2, ...
    }else{ out=out*wgt.ptb }
  }
  if(k0==1){out=out*yi}
  if(k0==2){out=c(out)*yi[,rep(1:py,py)]*yi[,rep(1:py,rep(py,py))]}
  as.matrix(sum.I(tt,"<=",xi,out))
}


PTB.Gamma.FUN.explicit <- function(wgt.ptb,data,gammahat, t0)
{
  logLamhat = gammahat[1:(length(t0))]; betahat = gammahat[-(1:length(t0))];
  xi = data[,1]; di = data[,2]; wgti = data[,3]; yi = data[,-(1:3),drop=F]; ebyi = c(exp(yi%*%betahat)); py = ncol(yi)
  tmpind = di==1; tj = xi[tmpind]; wgt.ptb.j = wgt.ptb[tmpind,]; pv = ncol(wgt.ptb)
  pi0.tj   = c(PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=0))/sum(wgti) 
  pi1.tj   =   PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=1)/sum(wgti)
  Ahat = PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=2)/sum(wgti)*pi0.tj - pi1.tj[,rep(1:py,py)]*pi1.tj[,rep(1:py,rep(py,py))]
  Ahat = matrix(apply(Ahat/pi0.tj^2,2,sum),ncol=py)
  
  term1 = t(wgt.ptb.j)%*%yi[tmpind,]
  term2 = t(wgt.ptb.j)%*%(pi1.tj/pi0.tj)
  
  pi1.tj.ptb = PI.k.FUN(tj,ebyi,xi,yi,wgt.ptb=wgt.ptb,k0=1)/sum(wgti) ## ny x nptb
  pi0.tj.ptb = PI.k.FUN(tj,ebyi,xi,yi,wgt.ptb=wgt.ptb,k0=0)/sum(wgti);## nptb x 1 
  term3 = t(rep(1,length(tj)))%*%((pi1.tj.ptb*pi0.tj-pi0.tj.ptb[,rep(1:pv,rep(py,pv))]*pi1.tj[,rep(1:py,pv)])/pi0.tj^2)
  term3 = t(matrix(term3,nrow=py))
  betastar = betahat + solve(Ahat)%*%t(term1 - term2 - term3)  # n.y x nptb
  
  term1 = sum.I(t0,">=",tj,wgt.ptb.j/pi0.tj)     # n.t0 x nptb
  term2 = sum.I(t0,">=",tj,pi0.tj.ptb/pi0.tj^2)  # n.t0 x nptb
  term3 = sum.I(t0,">=",tj,pi1.tj/pi0.tj^2)      # n.t0 x n.y
  term3 = term3%*%(betastar-betahat)
  logLamstar = logLamhat + (term1-term2-term3)/sum(wgti)/exp(logLamhat)
  rbind(logLamstar,betastar)
}



myfyk=function(yk, wgtk) {sum.I(yk, ">=", yk, wgtk)/sum(wgtk)}

Ptb.ROC.FUN <- function(ck,CondSck,Fck,uu0,type, yk)
{ 
  ooo = order(Fck)
  ck = ck[ooo]
  CondSck = CondSck[ooo]
  Fck = Fck[ooo]
  yk = yk[ooo]
  nc = length(ck)
  
  dFck = Fck - c(0,Fck[-nc])
  St0.Fck = cumsum(CondSck*dFck)## St0.Fck = P(T> t0,Y<=ck)
  Ft0.Fck = Fck-St0.Fck         ## Ft0.Fck = P(T<=t0,Y<=ck)
  St0 = max(St0.Fck)            ## St0     = P(T> t0      )
  Ft0 = 1-St0                   ## Ft0     = P(T<=t0      )
  FPR.c = (St0-St0.Fck)/St0     ## P(Y> ck|T> t0)
  TPR.c = (Ft0-Ft0.Fck)/Ft0     ## P(Y> ck|T<=t0)
  NPV.c = St0.Fck/Fck           ## P(T> t0|Y<=ck)
  PPV.c = (Ft0-Ft0.Fck)/(1-Fck) ## P(T<=t0|Y> ck)
  AUC = sum(TPR.c*(FPR.c-c(FPR.c[-1],0)))
  accuracy.out = data.frame("cutoff"=yk,"FPR"=FPR.c,"TPR"=TPR.c,"NPV"=NPV.c,"PPV"=PPV.c)
  junk = accuracy.out;
  
  ind0 = match(type,c("cutoff", "FPR","TPR","NPV","PPV")); 
  uuk = junk[,ind0]; junk = junk[order(uuk),-ind0]; uuk = sort(uuk); 
  tmpind = NULL; 
  for(i in 1:length(uu0)) {
    if(ind0==2){tpid = sum.I(uu0[i],">=",uuk)}else{tpid = sum.I(uu0[i],">",uuk)} 
    tmpind = c(tmpind,tpid)}
  
  
  c("AUC"=AUC,t(junk[tmpind,]))
}
