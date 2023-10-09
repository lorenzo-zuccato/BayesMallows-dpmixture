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
