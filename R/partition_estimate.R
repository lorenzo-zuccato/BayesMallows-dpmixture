partition_estimate <- function(model_fit){
  stopifnot(inherits(model_fit, "BayesMallowsDPMixture"))

  co_clus <- model_fit$co_clus

  if (is.null(co_clus)) {
    stop("model_fit$co_clus is missing. Please compute the co-clustering matrix usign function co_clustering.")
  }

  return(minVI(co_clus))
}
