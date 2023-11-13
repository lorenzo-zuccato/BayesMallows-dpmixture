co_clustering <- function(model_fit, burnin = model_fit$burnin){
  stopifnot(inherits(model_fit, "BayesMallowsDPMixture"))

  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  co_clus <- matrix(0, nrow = model_fit$n_assessors, ncol = model_fit$n_assessors)

  m <- model_fit$cluster_assignment
  for(iteration in m$iteration[m$iteration > burnin]){
    for(assessor1 in 2:fit$n_assessors){
      data_iter <- m[m$iteration == iteration , ]
      value <- data_iter$value[data_iter$assessor == assessor1]
      indices <- data_iter$assessor[data_iter$value == value]
      co_clus[assessor1, indices] <- co_clus[assessor1, indices] + 1
    }
  }

  return(co_clus / length(m$iteration[m$iteration > burnin]))
}

PROBLEMA:
  > max(fit$alpha$iteration)
[1] 4991
> max(fit$cluster_assignment$iteration)
[1] 500
> max(fit$rho$iteration)
[1] 4991
