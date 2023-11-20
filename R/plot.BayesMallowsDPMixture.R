#' Plot Alpha Posterior Distribution
#'
#' @param x An object of type \code{BayesMallowsDPMixture}, returned from
#'   \code{\link{compute_mallows_dpmixture}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{x$burnin}, and must be
#' provided if \code{x$burnin} does not exist. See \code{\link{assess_convergence_dpmixture}}.
#'
#' @export
plot.BayesMallowsDPMixture <- function(x, y = NULL, burnin = x$burnin, parameter = "n_clusters", ...) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  if (x$nmc <= burnin) stop("burnin must be <= nmc")

  if(parameter == "n_clusters"){
    df <- data.frame(n_clus = x$n_clusters[-seq(1, ceiling(burnin / x$clus_thin))])

    p <- ggplot2::ggplot(df, aes(n_clus)) +
              ggplot2::geom_bar(aes(y = (after_stat(count))/sum(after_stat(count))), color = "darkblue", fill = "darkblue")  +
              ggplot2::xlab("n_clusters") +
              ggplot2::ggtitle("Posterior distribution of number of occupied clusters") +
              ggplot2::scale_y_continuous(name = "Posterior probability",  breaks = c(0, 0.25, 0.5, 0.75, 1))

    return(p)
  }else if(parameter == "co_clustering"){
    co_clustering <- x$co_clustering

    if (is.null(co_clustering)) {
      stop("model_fit$co_clustering is missing. Please compute the co-clustering matrix usign function compute_co_clustering.")
    }

    assessor1 <- rep(seq(1, x$n_assessors), x$n_assessors)
    assessor2 <- rep(seq(1, x$n_assessors), each = x$n_assessors)
    probability <- c(co_clustering$co_clus_matrix)

    df <- data.frame(
      Assessor_1 = assessor1,
      Assessor_2 = assessor2,
      Probability = probability
    )

    p <- ggplot2::ggplot(data = df, aes(x = Assessor_1, y = Assessor_2)) +
      ggplot2::geom_raster(aes(fill = Probability)) +
      ggplot2::xlab("Assessors") +
      ggplot2::ylab("Assessors") +
      ggplot2::ggtitle("Co-clustering matrix") +
      ggplot2::scale_x_discrete(breaks = co_clustering$order) +
      ggplot2::scale_y_discrete(breaks = co_clustering$order)

    return(p)
  }else if(parameter == "alpha"){
    if (is.null(x$partition)) {
      stop("model_fit$partition is missing. Please compute the partition wrt which conditioning with function partition_estimate.")
    }

    df <- x$cluster_assignment
    df$partition <- rep(x$partition$cl, nrow(x$cluster_assignment) / x$n_assessors)
    df <- merge(x$alpha[(x$alpha$iteration > burnin) & is.finite(x$alpha$value), ],
                df[df$iteration > burnin , ],
                by.x = c("cluster", "iteration", "chain"), by.y = c("value", "iteration", "chain"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density(na.rm = TRUE) +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ylab("Posterior density")

    if (max(x$partition$cl) > 1) {
      p <- p + ggplot2::facet_wrap(~ .data$partition, scales = "free_x")
    }

    return(p)
  } else stop("parameter must be either \"n_clusters\", \"co_clustering\", or \"alpha\".")
}
