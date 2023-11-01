#include <RcppArmadillo.h>
#include "misc.h"
#include "mixtures.h"
#include "dpmixtures.h"
#include "distances.h"
#include "missing_data.h"
#include "pairwise_comparisons.h"
#include "parameterupdates.h"
#include "rmallows.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Worker function for computing the posterior distribution.
//'
//' @param rankings A set of complete rankings, with one sample per column.
//' With n_assessors samples and n_items items, rankings is n_items x n_assessors.
//' @param obs_freq  A vector of observation frequencies (weights) to apply to the rankings.
//' @param nmc Number of Monte Carlo samples.
//' @param constraints List of lists of lists, returned from `generate_constraints`.
//' @param cardinalities Used when metric equals \code{"footrule"} or
//' \code{"spearman"} for computing the partition function. Defaults to
//' \code{R_NilValue}.
//' @param logz_estimate Estimate of the log partition function.
//' @param metric The distance metric to use. One of \code{"spearman"},
//' \code{"footrule"}, \code{"kendall"}, \code{"cayley"}, or
//' \code{"hamming"}.
//' @param error_model Error model to use.
//' @param Lswap Swap parameter used by Swap proposal for proposing rank augmentations in the case of non-transitive pairwise comparisons.
//' @param leap_size Leap-and-shift step size.
//' @param alpha_prop_sd Standard deviation of proposal distribution for alpha.
//' @param alpha_init Initial value of alpha.
//' @param alpha_jump How many times should we sample \code{rho} between
//' each time we sample \code{alpha}. Setting \code{alpha_jump} to a high
//' number can significantly speed up computation time, since we then do not
//' have to do expensive computation of the partition function.
//' @param lambda Parameter of the prior distribution.
//' @param alpha_max Maximum value of \code{alpha}, used for truncating the exponential prior distribution.
//' @param psi Hyperparameter for the Dirichlet prior distribution used in clustering.
//' @param rho_thinning Thinning parameter. Keep only every \code{rho_thinning} rank
//' sample from the posterior distribution.
//' @param aug_thinning Integer specifying the thinning for data augmentation.
//' @param clus_thin Integer specifying the thinning for saving cluster assignments.
//' @param save_aug Whether or not to save the augmented data every
//' \code{aug_thinning}th iteration.
//' @param verbose Logical specifying whether to print out the progress of the
//' Metropolis-Hastings algorithm. If \code{TRUE}, a notification is printed every
//' 1000th iteration.
//' @param kappa_1 Hyperparameter for \eqn{theta} in the Bernoulli error model. Defaults to 1.0.
//' @param kappa_2 Hyperparameter for \eqn{theta} in the Bernoulli error model. Defaults to 1.0.
//'
// [[Rcpp::export]]
Rcpp::List run_mcmc_dpmixture(arma::mat rankings, arma::vec obs_freq, int nmc,
                            Rcpp::List constraints,
                            Rcpp::Nullable<arma::vec> cardinalities,
                            Rcpp::Nullable<arma::vec> logz_estimate,
                            Rcpp::Nullable<arma::mat> rho_init,
                            std::string metric = "footrule",
                            std::string error_model = "none",
                            int Lswap = 1,
                            int leap_size = 1,
                            double alpha_prop_sd = 0.5,
                            double alpha_init = 5,
                            int alpha_jump = 1,
                            double lambda = 0.1,
                            double alpha_max = 1e6,
                            int psi = 10,
                            int rho_thinning = 1,
                            int aug_thinning = 1,
                            int clus_thin = 1,
                            bool save_aug = false,
                            bool verbose = false,
                            double kappa_1 = 1.0,
                            double kappa_2 = 1.0){
  // The number of items ranked
  int n_items = rankings.n_rows;

  // The log factorial of the number of items is computed only once
  double log_fact_n_item = log_factorial(n_items);

  // The number of assessors
  int n_assessors = rankings.n_cols;

  bool augpair = (constraints.length() > 0);
  bool any_missing = !is_finite(rankings);

  umat missing_indicator;
  uvec assessor_missing;

  if(any_missing){
    // Converting to umat will convert NA to 0, but might cause clang-UBSAN error, so converting explicitly.
    rankings.replace(datum::nan, 0);
    missing_indicator = conv_to<umat>::from(rankings);
    missing_indicator.transform( [](int val) { return (val == 0) ? 1 : 0; } );
    assessor_missing = conv_to<uvec>::from(sum(missing_indicator, 0));
    initialize_missing_ranks(rankings, missing_indicator, assessor_missing);
  } else {
    missing_indicator.reset();
    assessor_missing.reset();
  }

  // Clustering
  int n_cluster_assignments = std::ceil(static_cast<double>(nmc * 1.0 / clus_thin));

  umat cluster_assignment(n_assessors, n_cluster_assignments);
  cluster_assignment.col(0) = initialize_cluster_assignment(n_assessors, psi);
  uvec current_cluster_assignment = cluster_assignment.col(0);

  uvec n_clusters(n_cluster_assignments);
  uvec current_clusters = unique(current_cluster_assignment);
  n_clusters(0) = current_clusters.n_elem;
  int current_n_clusters = n_clusters(0);
  unsigned int max_cluster_index = current_clusters.max();
  unsigned int min_cluster_index = 0;

  // Declare the cube to hold the latent ranks
  cube rho(n_items, current_n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
  rho.fill(datum::nan);
  rho.slice(0) = initialize_rho(n_items, current_n_clusters, rho_init);
  mat rho_old = rho(span::all, span::all, span(0));

  // Declare the vector to hold the scaling parameter alpha
  mat alpha(current_n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
  alpha.fill(datum::nan);
  alpha.col(0).fill(alpha_init);

  // If the user wants to save augmented data, we need a cube
  cube augmented_data;
  if(save_aug){
    augmented_data.set_size(n_items, n_assessors, std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
    augmented_data.slice(0) = rankings;
  }

  // Declare indicator vectors to hold acceptance or not
  mat alpha_acceptance(nmc - 1, current_n_clusters);
  alpha_acceptance.fill(datum::nan);
  mat rho_acceptance(nmc - 1, current_n_clusters);
  rho_acceptance.fill(datum::nan);

  vec aug_acceptance;
  if(any_missing | augpair){
    aug_acceptance = ones<vec>(n_assessors);
  } else {
    aug_acceptance.reset();
  }

  // Declare vector with Bernoulli parameter for the case of intransitive preferences
  vec theta, shape_1, shape_2;
  if(error_model == "bernoulli"){
    theta = zeros<vec>(nmc);
    shape_1 = zeros<vec>(nmc);
    shape_2 = zeros<vec>(nmc);
    shape_1(0) = kappa_1;
    shape_2(0) = kappa_2;
  } else {
    theta.reset();
    shape_1.reset();
    shape_2.reset();
  }

  // Other variables used
  int alpha_index = 0, rho_index = 0, aug_index = 0, cluster_assignment_index = 0;
  vec alpha_old = alpha.col(0);

  uvec element_indices = regspace<uvec>(0, rankings.n_rows - 1);

  // This is the Metropolis-Hastings loop

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0,
  // and this has been done above
  for(int t = 1; t < nmc; ++t){
    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << t << " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }
    }

    if(error_model == "bernoulli"){
      update_shape_bernoulli(shape_1(t), shape_2(t), kappa_1, kappa_2,
                             rankings, constraints);

      // Update the theta parameter for the error model, which is independent of cluster
      theta(t) = rtruncbeta(shape_1(t), shape_2(t), 0.5);
    }

    for(uvec::iterator i = current_clusters.begin(); i != current_clusters.end(); ++i){
      update_rho(rho, rho_acceptance, rho_old, rho_index, *i,
                 rho_thinning, alpha_old(*i), leap_size,
                 rankings.submat(element_indices, find(current_cluster_assignment == *i)),
                 metric, n_items, t, element_indices, obs_freq, min_cluster_index);
    }

    if(t % alpha_jump == 0) {
      ++alpha_index;
      for(uvec::iterator i = current_clusters.begin(); i != current_clusters.end(); ++i){
        alpha(*i, alpha_index) = update_alpha(alpha_acceptance, alpha_old(*i),
              rankings.submat(element_indices, find(current_cluster_assignment == *i)),
              obs_freq(find(current_cluster_assignment == *i)),
              *i, rho_old.col(*i), alpha_prop_sd, metric, lambda, t, cardinalities, logz_estimate, alpha_max);
      }
      // Update alpha_old
      alpha_old = alpha.col(alpha_index);
    }

    current_cluster_assignment = update_cluster_labels_dpmixture(rankings, rho, rho_old, rho_acceptance,
                                                                 alpha, alpha_old, alpha_acceptance, nmc,
                                                                 current_cluster_assignment, current_clusters, current_n_clusters,
                                                                 max_cluster_index, n_items, log_fact_n_item, lambda, alpha_max, psi,
                                                                 leap_size, metric, cardinalities, logz_estimate);

    if(t % clus_thin == 0){
      ++cluster_assignment_index;
      cluster_assignment.col(cluster_assignment_index) = current_cluster_assignment;
      n_clusters(cluster_assignment_index) = current_n_clusters;
    }

    // Perform data augmentation of missing ranks, if needed
    if(any_missing){
      update_missing_ranks(rankings, current_cluster_assignment, aug_acceptance, missing_indicator,
                             assessor_missing, alpha_old, rho_old, metric);
    }

    // Perform data augmentation of pairwise comparisons, if needed
    if(augpair){
      augment_pairwise(rankings, current_cluster_assignment, alpha_old, 0.1, rho_old,
                         metric, constraints, aug_acceptance, error_model, Lswap);
    }

    // Save augmented data if the user wants this. Uses the same index as rho.
    if(save_aug & (t % aug_thinning == 0)){
      ++aug_index;
      augmented_data.slice(aug_index) = rankings;
    }
  }

  vec alpha_acceptance_prob = ones(max_cluster_index);
  vec rho_acceptance_prob = ones(max_cluster_index);
  vec b;
  uvec temp_indices;

  for(unsigned int i = 0; i < max_cluster_index; ++i){
    b = alpha_acceptance.col(i);
    temp_indices = find_finite(b);
    alpha_acceptance_prob(i) = sum(b(temp_indices)) / temp_indices.n_elem;

    b = rho_acceptance.col(i);
    temp_indices = find_finite(b);
    rho_acceptance_prob(i) = sum(b(temp_indices)) / temp_indices.n_elem;
  }

  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance_prob,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance_prob,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("shape1") = shape_1,
    Rcpp::Named("shape2") = shape_2,
    Rcpp::Named("cluster_assignment") = cluster_assignment + 1,
    Rcpp::Named("n_clusters") = n_clusters,
    Rcpp::Named("augmented_data") = augmented_data,
    Rcpp::Named("any_missing") = any_missing,
    Rcpp::Named("augpair") = augpair,
    Rcpp::Named("aug_acceptance") = aug_acceptance / nmc,
    Rcpp::Named("n_assessors") = n_assessors,
    Rcpp::Named("obs_freq") = obs_freq
  );
}
