#include <RcppArmadillo.h>
#include "sample.h"
#include "partitionfuns.h"
#include "distances.h"
#include "misc.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

uvec initialize_cluster_assignment(int n_assessors, int psi){
    uvec cluster_assignment(n_assessors);
    cluster_assignment(0) = 0;
    for(int i = 1; i < n_assessors; ++i){
        vec prob = initial_assignment_prob(cluster_assignment.subvec(0,i-1), psi);

    }
}

vec initial_assignment_prob(const uvec already_assigned, int psi){
    int n_clusters = already_assigned.max()+1;
    vec prob(n_clusters+1);
    for(int i = 0; i < n_clusters; ++i){
        prob(i) = find(already_assigned == i).n_elem;
    }
    prob(n_clusters) = psi;
    return normalise(prob, 1);
}
