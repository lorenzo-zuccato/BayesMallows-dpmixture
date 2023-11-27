compute_co_clustering <- function(model_fit, burnin = model_fit$burnin){
  stopifnot(inherits(model_fit, "BayesMallowsDPMixture"))

  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  co_clus <- matrix(0, nrow = model_fit$n_assessors, ncol = model_fit$n_assessors)

  m <- model_fit$cluster_assignment[model_fit$cluster_assignment$iteration > burnin , ]
  for(iteration in unique(m$iteration)){
    for(assessor1 in 1:fit$n_assessors){
      data_iter <- m[m$iteration == iteration , ]
      value <- data_iter$value[data_iter$assessor == assessor1]
      indices <- data_iter$assessor[data_iter$value == value]
      co_clus[assessor1, indices] <- co_clus[assessor1, indices] + 1
    }
  }

  #h <- heatmap(co_clus)$rowInd
  #co_clus <- co_clus[h, h]

  #co_clustering <- list(co_clus_matrix = co_clus / length(unique(m$iteration)),
 #                       order = h)

  return(co_clus / length(unique(m$iteration)))
}
