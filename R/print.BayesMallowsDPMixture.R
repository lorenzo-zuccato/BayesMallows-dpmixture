#' Print Method for BayesMallowsDPMixture Objects
#'
#' The default print method for a \code{BayesMallowsDPMixture} object.
#'
#'
#' @param x An object of type \code{BayesMallowsDPMixture}, returned from
#'   \code{\link{compute_mallows_dpmixure}}.
#'
#' @param ... Other arguments passed to \code{print} (not used).
#'
#' @export
#'
#' @family posterior quantities
print.BayesMallowsDPMixture <- function(x, ...) {
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because print.BayesMallowsDPMixture must have the same
  # required arguments as base::print.

  if (is.null(x$n_items) || is.null(x$n_assessors)) {
    stop("BayesMallowsDPMixture object must have elements n_items and n_assessors.")
  }
  cat("Bayesian Mallows Model DPMixture with", x$n_items, "items and", x$n_assessors, "assessors.\n")
  cat("Use functions assess_convergence_dpmixture() or plot() to visualize the object.")
}
