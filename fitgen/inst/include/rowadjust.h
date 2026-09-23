#ifndef __ROWADJUST_H
#define __ROWADJUST_H

#include <fstream>

#include <Rcpp.h>
#include <diftimehelper.h>

#include <blissowen.h>
#include <blissowen_for_ext.h>
#include <stglm.h>

#define PHI_EST_NOT_NEEDED	0
#define PHI_EST_SINGLE		1
#define PHI_EST_COMMON		2
#define PHI_EST_GROUPED		3
#define PHI_EST_COMMON_EDGER	4
#define PHI_EST_TAGWISE_EDGER	5
#define N_PHI_EST		6

// Add constants and names (all the name in lower case) for more Phi estimation ways in
// future versions and increase the N_PHI_EST constant as needed
// ALSO, MODIFY the array procestnames in rowdajust.cpp
// It must be such that procestnames[CONSTANT] be the corresponding name

#endif

