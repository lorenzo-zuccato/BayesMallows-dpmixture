#include <RcppArmadillo.h>
#include "dpmixtures.h"
#include "rmallows.h"
#include "sample.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::uvec initialization(){
    uvec x(2);
    x(0) = 4000;
    uvec m(3);
    m(0) = 2;
    m(1) = 8;
    m(2) = 9;
    vec pr(3);
    pr(0) = 0.8;
    pr(1) = 0.1;
    pr(2) = 0.1;
    x(span(1)) = sample(m, 1, false, pr);
    return x;
}
