#include <RcppArmadillo.h>
#include "dpmixtures.h"
#include "rmallows.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat initialization(int n_samples){
    arma::mat init(4,n_samples);
    arma::vec rho0(4);
    rho0(0)=1;
    rho0(1)=2;
    rho0(2)=3;
    rho0(3)=4;
    init = rmallows(rho0, 1, n_samples, 10, 1);
    return init;
}
