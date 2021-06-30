#include <Rcpp.h>
#include <cmath> 
using namespace Rcpp;

// Function by https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/

// [[Rcpp::export]]
bool approximatelyEqualRel(double a, double b, double relEpsilon = 0.0001){
  double diff{ std::abs(a - b) };
  
  return (diff <= (std::max(std::abs(a), std::abs(b)) * relEpsilon));
}

// [[Rcpp::export]]
bool approximatelyEqualAbsRel(double a, double b, double absEpsilon = 0.0001, double relEpsilon = 0.0001){
  // Check if the numbers are really close -- needed when comparing numbers near zero.
  double diff{ std::abs(a - b) };
  if (diff <= absEpsilon)
    return true;
  
  // Otherwise fall back to Knuth's algorithm
  return (diff <= (std::max(std::abs(a), std::abs(b)) * relEpsilon));
}


/*** R
*/

