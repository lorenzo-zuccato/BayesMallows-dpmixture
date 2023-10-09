#ifndef DPMIXTURES_H
#define DPMIXTURES_H

#include "RcppArmadillo.h"

arma::vec initial_assignment_prob(const arma::uvec already_assigned, const int psi, const int n_clus);
arma::uvec initialize_cluster_assignment(const int n_assessors, const int psi);
#endif
