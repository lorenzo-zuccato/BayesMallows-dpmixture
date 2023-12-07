#' Trace Plots from Metropolis-Hastings Algorithm
#'
#' @param model_fit A fitted model object of class \code{BayesMallowsDPMixture} returned from
#'  \code{\link{compute_mallows_dpmixture}}.
#'
#' @param parameter Character string specifying which parameter to plot. Available
#' options are \code{"n_clusters"}, \code{"alpha"}, \code{"rho"}, \code{"empirical_cluster_probs"},
#'  \code{"Rtilde"} or \code{"theta"}.
#'
#' @param n Integer specifying the number of clusters to visualize. The first \code{n} clusters
#' ordered by number of iterations of persistence in the chain are chosen.
#'
#' @seealso \code{\link{compute_mallows_dpmixture}}
#'
#' @export
#' @family diagnostics
assess_convergence_dpmixture <- function(model_fit, parameter = "n_clusters", n = NULL, items = NULL, assessors = NULL) {
  stopifnot(inherits(model_fit, "BayesMallowsDPMixture"))

  if((parameter == "n_clusters") & !is.null(n)){
    stop("n must be specified only when parameter is set to \"alpha\", \"rho\" or \"empirical_cluster_probs\".")
  }

  if(parameter == "n_clusters"){
    m <- data.frame(n=model_fit$n_clusters)

    p <- ggplot(m, aes(x=unique(model_fit$cluster_assignment$iteration), y=n)) +
        geom_line(color="red")+
        ggplot2::xlab("Iteration") +
        ggplot2::ylab(expression(n)) +
        ggtitle("Number of non-empty clusters")

    return(p)

  } else if (parameter == "Rtilde") {

    trace_rtilde(model_fit, items, assessors)

  } else if (parameter == "theta") {

    trace_theta(model_fit)

  } else if ((parameter == "alpha") | (parameter == "empirical_cluster_probs") | (parameter == "rho")){
    if (is.null(n)) stop("Please specify the number of clusters to visualize.")
    if (n < 1) stop("The number of clusters to visualize must be at least 1.")

    persistence_clusters <- aggregate(iteration ~ value, data = model_fit$cluster_assignment, FUN=function(x) length(unique(x)))
    persistence_clusters <- persistence_clusters[order(-persistence_clusters$iteration), ]
    clusters <- persistence_clusters$value[1:n]

    if(parameter == "alpha"){
      m <- model_fit$alpha[model_fit$alpha$cluster %in% clusters ,]

      p <- ggplot2::ggplot(m, ggplot2::aes(
             x = iteration, y = value,
             group = cluster, color = cluster)) +
             ggplot2::xlab("Iteration") +
             ggplot2::ylab(expression(alpha)) +
             ggplot2::labs(color = "Cluster") +
             ggplot2::geom_line(ggplot2::aes(color = cluster)) +
             ggplot2::theme() +
             ggplot2::facet_wrap(ggplot2::vars(length(clusters)),
                            labeller = ggplot2::as_labeller(cluster_labeler_function),
                            scales = "free_y")
      return(p)

    }else if(parameter == "empirical_cluster_probs") {
      m <- model_fit$cluster_assignment[model_fit$cluster_assignment$value %in% clusters ,]

      m <- aggregate(m$chain, list(cluster = m$value, iteration =m$iteration), FUN = function(x){
                                                    return(length(x) / model_fit$n_assessors)})
      colnames(m) <- c("cluster", "iteration", "count")
      m2 <- model_fit$alpha[model_fit$alpha$cluster %in% clusters ,]

      m <- merge(m2, m, by=c("iteration", "cluster"), all = TRUE)
      m <- m[!is.na(m$count) , ]

      p <- ggplot2::ggplot(m, ggplot2::aes(
        x = iteration, y = count,
        group = cluster, color = cluster)) +
        ggplot2::xlab("Iteration") +
        ggplot2::ylab(expression(hat(tau))) +
        ggplot2::labs(color = "Cluster") +
        ggplot2::geom_line(ggplot2::aes(color = cluster)) +
        ggplot2::theme() +
        ggplot2::facet_wrap(ggplot2::vars(length(clusters)),
                            labeller = ggplot2::as_labeller(cluster_labeler_function),
                            scales = "free_y")
      return(p)

    } else if (parameter == "rho") {

      if (is.null(items) && model_fit$n_items > 5) {
        message("Items not provided by user. Picking 5 at random.")
        items <- sample.int(model_fit$n_items, 5)
      } else if (is.null(items) && model_fit$n_items > 0) {
        items <- seq.int(from = 1, to = model_fit$n_items)
      }

      if (!is.character(items)) {
        items <- model_fit$items[items]
      }

      df <- model_fit$rho[(model_fit$rho$item %in% items) & (model_fit$rho$cluster %in% clusters), , drop = FALSE]

      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$item)) +
        ggplot2::geom_line() +
        ggplot2::theme(legend.title = ggplot2::element_blank()) +
        ggplot2::xlab("Iteration") +
        ggplot2::ylab(expression(rho)) +
        ggplot2::facet_wrap(~factor(.data$cluster, clusters))

      return(p)

    } else stop("parameter must be either \"n_clusters\", \"alpha\", \"rho\", \"empirical_cluster_probs\", \"Rtilde\" or \"theta\".")
  }
}
