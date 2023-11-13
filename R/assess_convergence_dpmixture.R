#' Trace Plots from Metropolis-Hastings Algorithm
#'
#' @param model_fit A fitted model object of class \code{BayesMallowsDPMixture} returned from
#'  \code{\link{compute_mallows_dpmixture}}.
#'
#' @param parameter Character string specifying which parameter to plot. Available
#' options are \code{"n_clusters"}, \code{"alpha"} or \code{"empirical_cluster_probs"}.
#'
#' @param n Integer specifying the number of clusters to visualize. The first \code{n} clusters
#' ordered by number of iterations of persistence in the chain are chosen.
#'
#' @seealso \code{\link{compute_mallows_dpmixture}}
#'
#' @export
#' @family diagnostics
assess_convergence_dpmixture <- function(model_fit, parameter = "n_clusters", n = NULL) {
  stopifnot(inherits(model_fit, "BayesMallowsDPMixture"))

  if((parameter == "n_clusters") & !is.null(n)){
    stop(" n must be specified only when parameter is set to \"alpha\" or \"empirical_cluster_probs\".")
  }

  if(parameter == "n_clusters"){
    m <- model_fit$n_clusters

    plot(y = m, x = seq(1, model_fit$nmc, by = model_fit$clus_thin),
         main = "Number of non empty clusters", type = "l", col = "red",
         xlab = "Iteration", ylab = "")

  } else if((parameter == "alpha") | (parameter =="empirical_cluster_probs")){
    if (is.null(n)) stop("Please specify the number of clusters to visualize.")
    if (n < 1) stop("The number of clusters to visualize must be at least 1.")

    persistence_clusters <- fit$cluster_assignment %>%
      group_by(iteration, value) %>%
      filter(row_number() == 1)
    clusters <- names(sort(table(persistence_clusters$value), decreasing = TRUE)[1:n])

    if(parameter == "alpha"){
      m <- fit$alpha[fit$alpha$cluster %in% clusters ,]

      p <- ggplot2::ggplot(m, ggplot2::aes(
             x = iteration, y = value,
             group = cluster, color = cluster)) +
             ggplot2::xlab("Iteration") +
             ggplot2::ylab(expression(alpha)) +
             ggplot2::labs(color = "Cluster") +
             ggplot2::geom_line(ggplot2::aes(color = cluster)) +
             ggplot2::theme() +
             ggplot2::facet_wrap(ggplot2::vars(length(clusters)),
                            labeller = ggplot2::as_labeller(cluster_labeler_function), scales = "free_y")
    suppressWarnings(print(p))

    }else if(parameter == "empirical_cluster_probs") {
      m <- fit$cluster_assignment[fit$cluster_assignment$value %in% clusters ,]

      m <- aggregate(m, list(m$value, m$iteration), FUN = function(x){
                                                    return(length(x) / fit$n_assessors)})
      colnames(m) <- c("value", "iteration", "count")

      p <- ggplot2::ggplot(m, ggplot2::aes(
        x = iteration, y = count,
        group = value, color = value)) +
        ggplot2::xlab("Iteration") +
        ggplot2::ylab(expression(hat(tau))) +
        ggplot2::labs(color = "Value") +
        ggplot2::geom_line(ggplot2::aes(color = value)) +
        ggplot2::theme() +
        ggplot2::facet_wrap(ggplot2::vars(length(clusters)),
                            labeller = ggplot2::as_labeller(cluster_labeler_function), scales = "free_y")
      suppressWarnings(print(p))

    }
  }else stop("parameter must be either \"n_clusters\", \"alpha\" or \"empirical_cluster_probs\".")
}
