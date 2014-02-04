#ifndef _survAccuracyMeasuresRcpp_Code_h
#define _survAccuracyMeasuresRcpp_Code_h

#include <RcppArmadillo.h>

arma::uvec myRank( arma::colvec x);
arma::mat CSumI( arma::colvec yy, int FUN, arma::colvec Yi, arma::mat Vi, bool v);
arma::mat Vec2Mat(arma::colvec yy, int Nrows);
arma::colvec myPmin( arma::colvec myvec, double x) ;

RcppExport SEXP getWEXP(SEXP datar, SEXP Yr, SEXP N, SEXP RT_outr, SEXP predictTime, SEXP resid_sco, SEXP fitvar);

#endif
