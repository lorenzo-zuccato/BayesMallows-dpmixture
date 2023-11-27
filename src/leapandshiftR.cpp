#include <RcppArmadillo.h>
#include "leapandshift.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::vec leap_and_shiftR(const arma::vec& rho, int leap_size){
  arma::vec rho_proposal;
  arma::uvec indices;
  double prob_backward, prob_forward;

  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho, leap_size, true);

  return  rho_proposal;
}
