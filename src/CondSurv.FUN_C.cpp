// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "Code.h"

using namespace std; 
using namespace Rcpp;
using namespace arma; 

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
arma::mat myDnorm( arma::mat X){
  
  return as_scalar(1/sqrt(2.0*datum::pi))*exp(-pow(X, 2)/2.0); 
  
}

// [[Rcpp::export]]
arma::mat CondSurv_FUN_C(arma::mat IPW, arma::colvec xi, arma::uvec di, arma::colvec yi, double tt0, double bw){
  
  /*for(b in 1:B0) {
      
      Shat.yk.ptb[,b] = CondSurv.FUN(wgtk.ptb[,b],xk, dk, yk, t0, hhat.v); 
    
    }
  */
  /*
  nv = length(xi);
  
  #kerni.yy = Kern.FUN(c(yi),c(yi),bw)*IPW; ## nv x ny matrix
  kerni.yy <- (VTM(c(yi), length(yi))-c(yi))/bw
  kerni.yy <- dnorm(kerni.yy)/bw
  
  skern.yy = colSums(kerni.yy)
  tmpind = (xi<=tt0)&(di==1); tj = xi[tmpind]; nj = length(tj)
  pi.tj.yy = sum.I(tj,"<=",xi,kerni.yy)/VTM(skern.yy,nj) ## nj x ny matrix ##    
  dLam.tj.yy = kerni.yy[tmpind,]/pi.tj.yy/VTM(skern.yy,nj); 
  dLam.tj.yy[is.na(dLam.tj.yy)] = 0
  Shat.t0.yi = exp(-colSums(dLam.tj.yy))
  Shat.t0.yi    
  */
  
  
  int N = xi.n_elem; 
  int B0 = IPW.n_cols; 
  
  arma::mat output(N, B0); // called in Shat.yk.ptb in ncc_NPkernel.R
  arma::uvec tmpind = find((xi <= tt0)%(di==1));
   
  arma::colvec tj = xi.elem(tmpind);
  int nj = tj.n_elem;
  
  //define matrices out here first so I am not creating copies in the future
  arma::mat kerniyy(yi.n_rows, yi.n_rows); 
  arma::colvec skernyy(yi.n_rows); 
  arma::mat tmpDenom(nj, yi.n_rows); 
  arma::mat pitjyy(nj, yi.n_rows); 
  arma::mat dLamtjyy(nj, yi.n_rows); 
  arma::rowvec Shatt0yi(yi.n_rows); 
  
  for(int b=0; b < B0; b++){
 
     kerniyy = Vec2Mat(yi, yi.n_rows); 
     kerniyy.each_col() -= yi;
     kerniyy /= bw; 
     kerniyy = myDnorm(kerniyy)/bw;
     kerniyy.each_col() %= IPW.col(b);  //element wise multiplication 
     
     skernyy = sum(kerniyy, 0).t(); 
      
     tmpDenom = Vec2Mat(skernyy, nj); 
     pitjyy = CSumI(tj, 1, xi, kerniyy, TRUE)/tmpDenom; 
     dLamtjyy = (kerniyy.rows(tmpind)/pitjyy)/tmpDenom;
     Shatt0yi = exp(-sum(dLamtjyy, 0));
     
     output.col(b) = Shatt0yi.t(); 
   }
  
  return output; 
  
}

