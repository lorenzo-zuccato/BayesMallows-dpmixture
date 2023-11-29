#' @export
partition_estimate <- function(model_fit){
  stopifnot(inherits(model_fit, "BayesMallowsDPMixture"))

  co_clus <- model_fit$co_clustering

  if (is.null(co_clus)) {
    stop("model_fit$co_clustering is missing. Please compute the co-clustering matrix usign function co_clustering.")
  }

  partition <- minVI(co_clus)

  return(partition)
}
