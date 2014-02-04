// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "Code.h"

using namespace std; 
using namespace arma;
using namespace Rcpp; 

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
   return x.sort() * 2;
}

// [[Rcpp::export]]
NumericMatrix CondSurv_FUN_C(NumericVector IPW, NumericVector xi, IntegerVector di, NumericVector yi, double tt0, double bw){
  
  
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
  
  int nv = xi.size(); 
  mat kerni_yy = Vec2Mat(Rcpp::as<arma::vec >(yi), yi.size()); 
  
  
  return Rcpp::as<Rcpp::NumericMatrix>(wrap(kerni_yy)); 
  
}

