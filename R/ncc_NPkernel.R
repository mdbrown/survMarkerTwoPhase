#functions to calculate estimates and SE based on non parametric kernel smoothing


## data: xi, di, vi, yi; 
NCC.NP.kn.PTB.s.y <- function(data,V.IND,Iik0=Iik0.mat,wgtk.ptb=NULL,B0=500,Zmatch.ind,t0, m.match)
{     
  

  bw.power = 0.3
  xi = data[,1]; 
  di = data[,2]; 
  vi = data[,3]; 
  yi = data[,-(1:3),drop=F];
  yi=as.matrix(yi); 
  
  n.t0 = length(t0);
  NN = length(xi); 
  nv = sum(vi);
  
  tj = V.IND[,1]; 
  pi.tj = V.IND[,2]/NN
  wgtk = 1/FNCC.WGT.FUN(data,V.IND,Iik0=Iik0.mat,Zmatch.ind, m.match = m.match); 
  
  xk = xi[vi==1]; dk = di[vi==1]; yk = yi[vi==1,,drop=F]; idk = (1:NN)[vi==1];  
  yk.old <- yk
  #wgtk.ptb = PtbNCC.WGT.FUN(data,V.IND,Iik0,Zmatch.ind=Zmatch.ind,B0=B0); wgtk.ptb = wgtk.ptb[match(idk,wgtk.ptb[,1]),-1]
  wgtk.ptb = PtbNCC.WGT.FUN(data,V.IND,Iik0,B0=B0, m.match=m.match); 
  wgtk.ptb = wgtk.ptb[match(idk,wgtk.ptb[,1]),-1]
  npy = dim(yk)[2]

  #if (npy>1) {
  #  betahat = coxph(Surv(xk,dk)~yk,weights=wgtk,robust=T); 
  #  betahat=betahat$coef; 
  #  ebyk = c(exp(yk%*%betahat));
    
  #  betaptb = PTB.beta.FUN.explicit(wgtk.ptb,cbind(xk,dk,wgtk,yk),betahat=betahat) 
  #  yk.ptb= matrix(0,nv,B0)
  #  yk.ptb = yk%*%betaptb
  #  yk = ebyk   	
  #}
  
  hhat.v = 1.06*min(sd(yk),IQR(yk)/1.34)*nv^(-bw.power)
  Shat.yk = CondSurv.FUN(wgtk, xk, dk, yk, t0, hhat.v); 

  Shat.yk.ptb <- CondSurv_FUN_C(wgtk.ptb,as.matrix(xk), dk, yk, t0, hhat.v)

  #if (npy>1) {
  #  for (b in 1:B0) {
  #    hhat.v = 1.06*min(sd(yk.ptb[,b]),IQR(yk.ptb[,b])/1.34)*nv^(-bw.power)
  #    Shat.yk.ptb[,b] = CondSurv.FUN(wgtk.ptb[,b],xk, dk, yk.ptb[,b], t0, hhat.v); }
  #} else {
  
  
 # Shat.yk.ptb2=matrix(0,nv,B0)

  #  for (b in 1:B0) {
  #Shat.yk.ptb2[,b] =  CondSurv.FUN(wgtk.ptb[,b],xk, dk, yk, t0, hhat.v); 
  #  }
  
  #}  
  
  Fyk = sum.I(yk, ">=", yk, wgtk)/sum(wgtk) 
  Fyk.ptb = matrix(0,length(yk),B0)
  
  for (b in 1:B0) {Fyk.ptb[,b] =sum.I(yk, ">=", yk, wgtk.ptb[,b])/sum(wgtk.ptb[,b])}
  #Fyk.ptb = apply(wgtk.ptb,2,myfyk(yk,wgtk))
  if (npy>1) {
    list(ck =yk, 
         beta.est=betahat,
         beta.sd =apply(betaptb,1,sd), 
         CondSck = Shat.yk,
         CondSck.ptb = Shat.yk.ptb, 
         Fck = Fyk,
         Fck.ptb = Fyk.ptb, 
         yk = yk.old)
  } else{ 
    list(ck =yk, 
         CondSck = Shat.yk,
         CondSck.ptb = Shat.yk.ptb, 
         Fck = Fyk,
         Fck.ptb = Fyk.ptb, 
         yk = yk.old)  }
}


PTB.beta.FUN.explicit <- function(wgt.ptb,data,betahat)
{ 
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
  betastar
}

CondSurv.FUN<-function(IPW, xi, di, yi, tt0, bw)
{
  nv = length(xi); 
  
  #kerni.yy = Kern.FUN(c(yi),c(yi),bw)*IPW; ## nv x ny matrix
 # tmp <- CondSurv_FUN_C(as.matrix(IPW), xi, di, yi, tt0, bw)

  kerni.yy <- (VTM(c(yi), length(yi))-c(yi))/bw
  kerni.yy <- (dnorm(kerni.yy)/bw)*IPW
  
  skern.yy = colSums(kerni.yy)
  tmpind = (xi<=tt0)&(di==1); tj = xi[tmpind]; nj = length(tj)
  pi.tj.yy = sum.I(tj,"<=",xi,kerni.yy)/VTM(skern.yy,nj) ## nj x ny matrix ##    
  dLam.tj.yy = kerni.yy[tmpind,]/pi.tj.yy/VTM(skern.yy,nj); 
  dLam.tj.yy[is.na(dLam.tj.yy)] = 0
  Shat.t0.yi = exp(-colSums(dLam.tj.yy))
  Shat.t0.yi    
}


Kern.FUN <- function(zz,zi,bw,kern0="gauss") ## returns an (n x nz) matrix ##
{ 
  
  out = (VTM(zz,length(zi))- zi)/bw
  switch(kern0,
         "epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
         "gauss"= dnorm(out)/bw)
}
