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
    const uvec& current_cluster_assignment;
    const uvec& current_clusters;
    const int& current_n_clusters;
    const vec& alpha_old,
    const int& n_items,
    const int& t,
    const std::string& metric,
    const Rcpp::Nullable<vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<vec> logz_estimate = R_NilValue,
){
  int n_assessors = dist_mat.n_rows;
  long n_item_factorial = factorial(n_items);
  uvec new_cluster_assignment(n_assessors);
  int n_in_cluster, cluster_index;
  vec assignment_probabilities;

  for(int i = 0; i < n_assessors; ++i){
    assignment_probabilities.set_size(current_n_clusters + 1);
    for(int j = 0; j < current_n_clusters; ++j){
        cluster_index = current_clusters(i);
        n_in_cluster = n_elem_in_cluster(current_cluster_assignment, n_assessors, cluster_index, i);
        assignment_probabilities(i) = std::log(n_in_cluster) - std::log(psi + n_assessors - 1) -
          get_partition_function(n_items, alpha_old(cluster_index), cardinalities, logz_estimate, metric) -
          alpha.old(cluster_index) / n_items * //!SOME DISTANCE HERE MISSING;
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

int n_elem_in_cluster(const uvec& current_cluster_assignment,
                      const int& n_assessors,
                      const int& cluster_index,
                      const int& i){
  int c = 0;
  for(int j = 0; j < n_assessors; ++j){
    if(current_cluster_assignment(j) == cluster_index & j != i) c++;
  }
  return c;
}

long factorial(const int n){
    long f = 1;
    for (int i=1; i<=n; ++i)
        f *= i;
    return f;
}
