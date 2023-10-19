#include <RcppArmadillo.h>
#include "sample.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


vec initial_assignment_prob(const uvec already_assigned, const int psi, const int n_clus){
    vec prob(n_clus + 1);
    for(int i = 0; i < n_clus; ++i){
        uvec elem_i = find(already_assigned == i);
        prob(i) = elem_i.n_elem;
    }
    prob(n_clus) = psi;
    return normalise(prob, 1);
}

uvec initialize_cluster_assignment(const int n_assessors, const int psi){
    uvec cluster_assignment(n_assessors);
    cluster_assignment(0) = 0;
    int n_clus = 1;
    for(int i = 1; i < n_assessors; ++i){
        n_clus = cluster_assignment.subvec(0, i - 1).max() + 1;
        vec assignment_prob(n_clus+1);
        assignment_prob = initial_assignment_prob(cluster_assignment.subvec(0, i - 1), psi, n_clus);
        cluster_assignment(span(i)) = sample(regspace<uvec>(0, n_clus), 1, false, assignment_prob);
    }
    return cluster_assignment;
}

uvec update_cluster_labels_dpmixture(
    //const mat& dist_mat,
    cube& rho,
    mat& rho_old,
    mat& alpha,
    vec& alpha_old,
    const uvec& current_cluster_assignment,
    uvec& current_clusters,
    int& current_n_clusters,
    int& max_cluster_index,
    const int& n_items,
    const double& log_fact_n_items,
    const int& t,
    const int& psi,
    const std::string& metric,
    const Rcpp::Nullable<vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<vec> logz_estimate = R_NilValue,
){
  int n_assessors = dist_mat.n_rows, cluster_index;
  uvec new_cluster_assignment(n_assessors);
  uvec n_in_cluster, possible_clusters;
  vec assignment_probabilities, part_fun, alpha;

  for(int i = 0; i < n_assessors; ++i){
    //Resizing objects based on current number of clusters
    assignment_probabilities.set_size(current_n_clusters + 1);
    n_in_cluster.set_size(current_n_clusters + 1);
    possible_clusters.set_size(current_n_clusters + 1);
    part_fun.set_size(current_n_clusters);
    alpha.set_size(current_n_clusters)

    // Calculating assignment probabilities
    n_in_cluster.subvec(0, current_n_clusters - 1) = cluster_count_excluding_i(current_cluster_assignment,cluster_index,
                                                                                n_assessors, i);
    for(int j = 0; j < current_n_clusters; ++j){
        cluster_index = current_clusters(j);
        part_fun(j) = get_partition_function(n_items, alpha_old(cluster_index), cardinalities, logz_estimate, metric);
        alpha(j) = alpha_old(cluster_index);
        n_in_cluster(j) = cluster_count_excluding_i(current_cluster_assignment,cluster_index,
                                                                                n_assessors, i);
    }
    // Compute the logarithm of the unnormalized probability
    assignment_probabilities.subvec(0, current_n_clusters - 1) = std::log(n_in_cluster) - std::log(psi + n_assessors -1)
                                        - part_fun - alpha_old_temp / n_items * DISTANCE//!SOME DISTANCE HERE MISSING;
    assignment_probabilities(current_n_clusters) = std::log(psi) - log_fact_n_items - std::log(psi + n_assessors - 1);
    // Exponentiate to get unnormalized prob relative to max
    vec probs = exp(assignment_probabilities - max(assignment_probabilities));
    assignment_probabilities = normalise(prob, 1);

    // Setting up possible cluster labels
    possible_clusters.subvec(0, current_n_clusters - 1) = current_clusters;
    possible_clusters(current_n_clusters) = max_cluster_index + 1;

    // Sampling new cluster for assessor i
    new_cluster_assignment(i) = sample(possible_clusters, 1, false, assignment_probabilities);

    // Updating current clusters
    current_clusters = unique(new_cluster_assignment);
    current_n_clusters = current_clusters.n_elem;

    if(new_cluster_assignment(i) == max_cluster_index + 1){
        max_cluster_index++;

        rho.resize(n_items, max_cluster_index + 1, rho.n_slices);
        rho(span::all, span(max_cluster_index ), span::all).fill(datum::nan);
        rho_old.resize(n_items, max_cluster_index + 1);

        alpha.resize(max_cluster_index - 1, alpha.n_cols);
        alpha.row(max_cluster_index).fill(datum::nan);
        alpha_old.resize(max_cluster_index + 1);

        alpha_old(max_cluster_index) =
        rho_old.col(max_cluster_index) = rmallows()
    }
  }

 //!DA QUA NON C'ENTRA. HO SOLO COPIATO LA FUNZIONE PER LE MISTURE FINITE

  for(int i = 0; i < n_clusters; ++i){
    // Compute the logarithm of the unnormalized probability
    assignment_prob.col(i) = std::log(cluster_probs(i)) -
      alpha_old(i) / n_items * dist_mat.col(i) -
      get_partition_function(n_items, alpha_old(i), cardinalities, logz_estimate, metric);
  }

  for(int i = 0; i < n_assessors; ++i){
    // Exponentiate to get unnormalized prob relative to max
    rowvec probs = exp(assignment_prob.row(i) -
      max(assignment_prob.row(i)));

    // Normalize with 1-norm
    assignment_prob.row(i) = normalise(probs, 1);
    new_cluster_assignment(span(i)) = sample(regspace<uvec>(0, probs.n_elem - 1), 1, false, assignment_prob.row(i).t());
  }

  if(save_ind_clus){
    assignment_prob.save(std::string("cluster_probs") + std::to_string(t + 1) + std::string(".csv"), csv_ascii);
  }
  return(new_cluster_assignment);
}

int cluster_count_excluding_i(const uvec& current_cluster_assignment,
                      const int& cluster_index,
                      const int& n_assessors,
                      const int& i){
  int c = 0;
  for(int j = 0; j < n_assessors; ++j){
            if(current_cluster_assignment(j) == cluster_index & j != i) c++;
        }
  return c;
}



double  log_factorial(const int n){
    double f = 0;
    for (int i=1; i<=n; ++i)
        f += std::log(i);
    return f;
}
