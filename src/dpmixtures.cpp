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
    int n_clus = 1;
    for(int i = 1; i < n_assessors; ++i){
        n_clus = cluster_assignment.subvec(0, i - 1).max() + 1;
        vec assignment_prob = initial_assignment_prob(cluster_assignment.subvec(0, i - 1), psi, n_clus);
        cluster_assignment(i) = sample(regspace<uvec>(0, n.clus), 1, false, assignment_prob);
    }
    return cluster_assignment;
}

vec initial_assignment_prob(const uvec already_assigned, int psi, int n_clus){
    vec prob(n_clus + 1);
    //for(int i = 0; i < n_clusters; ++i){
    //    prob(i) = find(already_assigned == i).n_elem;
    //}
    //prob(n_clus) = psi;
    //return normalise(prob, 1);
    for(int i = 0; i < n_clusters; ++i){
        prob(i) = find(already_assigned == i).n_elem / (i + psi);
    }
    prob(n_clus) = psi / (i + psi);
    return prob;
}
