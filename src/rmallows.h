#ifndef RMALLOWS_H
#define RMALLOWS_H

#include "RcppArmadillo.h"
arma::mat rmallows(
    arma::vec rho0,
    arma::vec obs_freq,
    double alpha0,
    int n_samples,
    int burnin,
    int thinning,
    int leap_size = 1,
    std::string metric = "footrule"
  );

#endif
