#include <RcppArmadillo.h>
#include "dpmixtures.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::uvec initialization(int n, int psi){
    arma::uvec init(n);
    init = initialize_cluster_assignment(n, psi);
    return init;
}
