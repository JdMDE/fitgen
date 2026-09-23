#ifndef __STGLM_H
#define __STGLM_H

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <algorithm>
#include <cctype>

#include "debug.h"

// Constants to identify the distributions
#define DPOISS     0 // Poisson
#define DBINLOGIT  1 // Binomial with logistic link function (ungrouped)
#define DNEGBIN    2 // Negative binomial
#define DBINLOGITW 3 // Grouped binomial with logistic link function 
#define NPOSSDIST  4 // Number of possible distributions

// Add constants and names (all the name in lower case) for more distributions in
// future versions and increase the NPOSSDIST constant as needed.
// ALSO, MODIFY the array distnames in stglm.cpp.
// It must be such that distnames[CONSTANT] be the corresponding name

int CheckDistName(std::string model);

Rcpp::List StGLM(Rcpp::NumericVector y,Rcpp::NumericMatrix X,std::string model,Rcpp::NumericVector offset,bool warn_me,double Phi,bool retmu);

#endif

