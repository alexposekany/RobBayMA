#ifndef APPROX_EQUAL_CPP_H
#define APPROX_EQUAL_CPP_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
bool approximatelyEqualRel(double a, double b, double relEpsilon = 0.0001);
  
#endif
  