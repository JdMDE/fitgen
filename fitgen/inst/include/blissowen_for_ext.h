#ifndef __BLISSOWEN_FOR_EXT_H
#define __BLISSOWEN_FOR_EXT_H

#include <Rcpp.h>

Rcpp::List BlissOwen(Rcpp::NumericVector y);
Rcpp::List CommonBlissOwen(Rcpp::NumericMatrix y,int maxit=100,double tol=1e-9,double posthres=0.0);
Rcpp::List GroupedCommonBlissOwen(Rcpp::NumericMatrix y,int maxit=100,double tol=1e-9,double tolphi=-1.0,double posthres=0.0,int down_ndev=3,int up_ndev=2,double down_maxprop=1.0,double up_maxprop=1.0,int niter=2);

#endif

