#' Plot Alpha Posterior Distribution
#'
#' @param x An object of type \code{BayesMallowsDPMixture}, returned from
#'   \code{\link{compute_mallows_dpmixture}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{x$burnin}, and must be
#' provided if \code{x$burnin} does not exist. See \code{\link{assess_convergence_dpmixture}}.
#'
#' @param parameter Character string defining the parameter to plot. Available
#' options are \code{"n_clusters"}, \code{"co_cluseting"}, \code{"alpha"}, \code{"rho"},
#' and \code{"theta"}.
#'
#' @param items The items to study in the diagnostic plot for \code{rho}. Either
#'   a vector of item names, corresponding to \code{x$items} or a
#'   vector of indices. If NULL, five items are selected randomly.
#'   Only used when \code{parameter = "rho"}.
#'
#' @param ... Other arguments passed to \code{plot} (not used).
#'
#' @export
plot.BayesMallowsDPMixture <- function(x, burnin = x$burnin, parameter = "n_clusters", items = NULL, ...) {
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
  } else if(parameter == "co_clustering"){
    co_clustering <- x$co_clustering

    if (is.null(co_clustering)) {
      stop("model_fit$co_clustering is missing. Please compute the co-clustering matrix usign function compute_co_clustering.")
    }

    assessor1 <- rep(seq(1, x$n_assessors), x$n_assessors)
    assessor2 <- rep(seq(1, x$n_assessors), each = x$n_assessors)

    order <- hclust(dist(co_clustering))$order

    probability <- c(co_clustering[order, order])

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
      ggplot2::scale_x_discrete(breaks = order) +
      ggplot2::scale_y_discrete(breaks = order)

    return(p)
  } else if(parameter == "alpha"){
    if (is.null(x$partition)) {
      stop("model_fit$partition is missing. Please compute the partition wrt which conditioning with function partition_estimate.")
    }

    df <- x$cluster_assignment
    df$partition <- rep(paste0("Cluster ", x$partition$cl), nrow(x$cluster_assignment) / x$n_assessors)
    df <- merge(x$alpha[(x$alpha$iteration > burnin), ],
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
  } else if (parameter == "rho"){
    if (is.null(x$partition)) {
      stop("model_fit$partition is missing. Please compute the partition wrt which conditioning with function partition_estimate.")
    }

    if (is.null(items) && x$n_items > 5) {
      message("Items not provided by user. Picking 5 at random.")
      items <- sample.int(x$n_items, 5)
    } else if (is.null(items) && x$n_items > 0) {
      items <- seq.int(from = 1, to = x$n_items)
    }

    if (!is.character(items)) {
      items <- x$items[items]
    }

    df <- x$cluster_assignment
    df$partition <- rep(paste0("Cluster ", x$partition$cl), nrow(x$cluster_assignment) / x$n_assessors)
    df <- merge(x$rho[x$rho$iteration > burnin & x$rho$item %in% items, , drop = FALSE],
                df[df$iteration > burnin , ],
                by.x = c("cluster", "iteration", "chain"), by.y = c("value", "iteration", "chain"))

    df <- aggregate(
      list(n = df$iteration),
      list(cluster = df$partition, item = df$item, value = df$value),
      FUN = length
    )

    df$pct <- ave(df$n, df$cluster, df$item, FUN = function(x){x / sum(x)})

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$pct)) +
      ggplot2::geom_col() +
      ggplot2::scale_x_continuous(labels = scalefun) +
      ggplot2::xlab("rank") +
      ggplot2::ylab("Posterior probability")

    if (max(x$partition$cl) == 1) {
      p <- p + ggplot2::facet_wrap(~ .data$item)
    } else {
      p <- p + ggplot2::facet_wrap(~ .data$cluster + .data$item)
    }

    return(p)
  } else if(parameter == "theta"){
    if (is.null(x$theta)) {
      stop("Please run compute_mallows with error_model = 'bernoulli'.")
    }

    df <- x$theta[x$theta$iteration > burnin, , drop = FALSE]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(theta)) +
      ggplot2::ylab("Posterior density")

    return(p)
  } else stop("parameter must be either \"n_clusters\", \"co_clustering\", \"alpha\", \"rho\" or \"theta\".")
}
