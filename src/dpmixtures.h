#ifndef DPMIXTURES_H
#define DPMIXTURES_H

#include "RcppArmadillo.h"

arma::vec initial_assignment_prob(const arma::uvec already_assigned, const int psi, const int n_clus);
arma::uvec initialize_cluster_assignment(const int n_assessors, const int psi);

int cluster_count_excluding_i(const arma::uvec& current_cluster_assignment,
                              const unsigned int& cluster_index,
                              const int& n_assessors,
                              const int& i);
double  log_factorial(const int n);

arma::uvec update_cluster_labels_dpmixture(const arma::mat& rankings,
                                    arma::cube& rho,
                                    arma::mat& rho_old,
                                    arma::mat& rho_acceptance,
                                    arma::mat& alpha,
                                    arma::vec& alpha_old,
                                    arma::mat& alpha_acceptance,
                                    const int& nmc,
                                    const arma::uvec& current_cluster_assignment,
                                    arma::uvec& current_clusters,
                                    int& current_n_clusters,
                                    unsigned int& max_cluster_index,
                                    const int& n_items,
                                    const double& log_fact_n_items,
                                    const double& lambda,
                                    const double& alpha_max,
                                    const int& psi,
                                    const int& leap_size,
                                    const std::string& metric,
                                    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                                    const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue);
#endif
