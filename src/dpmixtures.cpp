#include <RcppArmadillo.h>
#include "sample.h"
#include "partitionfuns.h"
#include "distances.h"
#include "misc.h"

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

int cluster_count_excluding_i(const uvec& current_cluster_assignment,
                      const unsigned int& cluster_index,
                      const int& n_assessors,
                      const int& i){
  int c = 0;
  for(int j = 0; j < n_assessors; ++j){
            if((current_cluster_assignment(j) == cluster_index) & (j != i)) c++;
        }
  return c;
}

double  log_factorial(const int n){
    double f = 0;
    for (int i=1; i<=n; ++i)
        f += std::log(i);
    return f;
}

uvec update_cluster_labels_dpmixture(const mat& rankings
                                    const vec& obs_freq,
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
                                    const double& lambda,
                                    const double& alpha_max,
                                    const int& psi,
                                    const int& leap_size,
                                    const std::string& metric,
                                    const Rcpp::Nullable<vec> cardinalities = R_NilValue,
                                    const Rcpp::Nullable<vec> logz_estimate = R_NilValue){
  int n_assessors = dist_mat.n_rows;
  bool disappearing_cluster = false;
  unsigned int cluster_index, disappearing_cluster_index;
  uvec new_cluster_assignment(n_assessors);
  uvec n_in_cluster, possible_clusters;
  vec assignment_probabilities, part_fun, alpha_old_temp, dist;
  mat rho_old_temp;

  for(int i = 0; i < n_assessors; ++i){
    //Resizing objects based on current number of clusters
    assignment_probabilities.set_size(current_n_clusters + 1);
    n_in_cluster.set_size(current_n_clusters + 1);
    possible_clusters.set_size(current_n_clusters + 1);
    part_fun.set_size(current_n_clusters);
    alpha_old_temp.set_size(current_n_clusters);
    rho_old_temp.set_size(n_items, current_n_clusters);

    // Calculating assignment probabilities
    for(int j = 0; j < current_n_clusters; ++j){
        cluster_index = current_clusters(j);
        part_fun(j) = get_partition_function(n_items, alpha_old(cluster_index), cardinalities, logz_estimate, metric);
        n_in_cluster(j) = cluster_count_excluding_i(current_cluster_assignment, cluster_index, n_assessors, i);
        // Checking if the cluster that i currently belongs to is going to disappear in this iteration
        if((cluster_index == current_cluster_assignment(i)) & (n_in_cluster(j) == 0)){
            disappearing_cluster = true;
            disappearing_cluster_index = cluster_index;
        }
        alpha_old_temp(j) = alpha_old(cluster_index);
        rho_old_temp.col(j) = rho_old.col(cluster_index);
    }
    // Calculating distances of assessor i from current cluster consensus
    dist = rank_dist_vec(rho_old_temp, rankings.col(i), metric, obs_freq); //!CHECK THAT THIS FUNCTION IS CORRECT, STRANGE BEHAVIOUR WITH OBS_FREQ
    // Compute the logarithm of the unnormalized probability
    assignment_probabilities.subvec(0, current_n_clusters - 1) = std::log(n_in_cluster) - std::log(psi + n_assessors -1)
                                        - part_fun - alpha_old_temp / n_items * dist;
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

    // Resize and sampling new parameters when a new cluster is created
    if(new_cluster_assignment(i) == max_cluster_index + 1){
        max_cluster_index++;

        rho.resize(n_items, max_cluster_index + 1, rho.n_slices);
        rho(span::all, span(max_cluster_index ), span::all).fill(datum::nan);
        rho_old.resize(n_items, max_cluster_index + 1);

        alpha.resize(max_cluster_index - 1, alpha.n_cols);
        alpha.row(max_cluster_index).fill(datum::nan);
        alpha_old.resize(max_cluster_index + 1);

        alpha_old(max_cluster_index) = rtruncexp(lambda, alpha_max);
        rho_old.col(max_cluster_index) = rmallows(rankings.col(i), obs_freq, alpha_old(max_cluster_index),
                                                  1, 1000*std::log(n_items), 0, leap_size, metric);
    }

    // Deleting parameters of the disappearing cluster
    if(disappearing_cluster){
        rho_old.col(disappearing_cluster).fill(datum::nan);
        alpha_old(disappearing_cluster) = datum::nan;
    }
  }
  return new_cluster_assignment;
}
