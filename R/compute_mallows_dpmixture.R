compute_mallows_dpmixture <- function(rankings = NULL,
                            preferences = NULL,
                            obs_freq = NULL,
                            metric = "footrule",
                            error_model = NULL,
                            clus_thin = 1L,
                            nmc = 2000L,
                            leap_size = max(1L, floor(n_items / 5)),
                            swap_leap = 1L,
                            rho_init = NULL,
                            rho_thinning = 1L,
                            alpha_prop_sd = 0.1,
                            alpha_init = 1,
                            alpha_jump = 1L,
                            lambda = 0.1,
                            alpha_max = 100,
                            psi = 0.1,
                            psi_init = 5,
                            save_aug = FALSE,
                            aug_thinning = 1L,
                            logz_estimate = NULL,
                            verbose = FALSE,
                            validate_rankings = TRUE,
                            na_action = "augment",
                            constraints = NULL,
                            save_ind_clus = FALSE,
                            seed = NULL
                            ) {
  if (!is.null(seed)) set.seed(seed)

  # Check if there are NAs in rankings, if it is provided
  if (!is.null(rankings)) {
    if (na_action == "fail" && any(is.na(rankings))) {
      stop("rankings matrix contains NA values")
    }

    if (na_action == "omit" && any(is.na(rankings))) {
      keeps <- apply(rankings, 1, function(x) !any(is.na(x)))
      print(paste("Omitting", sum(keeps), "rows from rankings due to NA values"))
      rankings <- rankings[keeps, , drop = FALSE]
    }
  }

  # Check that at most one of rankings and preferences is set
  if (is.null(rankings) && is.null(preferences)) {
    stop("Either rankings or preferences (or both) must be provided.")
  }

  if (is.null(preferences) && !is.null(error_model)) {
    stop("Error model requires preferences to be set.")
  }

  # Check if obs_freq are provided
  if (!is.null(obs_freq)) {
    if (is.null(rankings)) {
      stop("rankings matrix must be provided when obs_freq are provided")
    }
    if (nrow(rankings) != length(obs_freq)) {
      stop("obs_freq must be of same length as the number of rows in rankings")
    }
  }

  if (!swap_leap > 0) stop("swap_leap must be strictly positive")
  if (nmc <= 0) stop("nmc must be strictly positive")

  # Check that we do not jump over all alphas
  if (alpha_jump >= nmc) stop("alpha_jump must be strictly smaller than nmc")

  # Check that we do not jump over all rhos
  if (rho_thinning >= nmc) stop("rho_thinning must be strictly smaller than nmc")
  if (aug_thinning >= nmc) stop("aug_thinning must be strictly smaller than nmc")

  if (lambda <= 0) stop("exponential rate parameter lambda must be strictly positive")

  # Check that all rows of rankings are proper permutations
  if (!is.null(rankings) && validate_rankings && !all(apply(rankings, 1, validate_permutation))) {
    stop("invalid permutations provided in rankings matrix")
  }


  # Deal with pairwise comparisons. Generate rankings compatible with them.
  if (!is.null(preferences) && is.null(error_model)) {
    if (!inherits(preferences, "BayesMallowsTC")) {
      message("Generating transitive closure of preferences.")
      # Make sure the preference columns are double
      preferences$bottom_item <- as.numeric(preferences$bottom_item)
      preferences$top_item <- as.numeric(preferences$top_item)
      preferences <- generate_transitive_closure(preferences)
    }
    if (is.null(rankings)) {
      message("Generating initial ranking.")
      rankings <- generate_initial_ranking(preferences)
    }
  } else if (!is.null(error_model)) {
    stopifnot(error_model == "bernoulli")
    n_items <- max(c(preferences$bottom_item, preferences$top_item))
    n_assessors <- length(unique(preferences$assessor))
    if (is.null(rankings)) {
      rankings <- replicate(n_assessors, sample(x = n_items, size = n_items), simplify = "numeric")
      rankings <- matrix(rankings, ncol = n_items, nrow = n_assessors, byrow = TRUE)
    }
  }

  # Find the number of items
  n_items <- ncol(rankings)

  # If any row of rankings has only one missing value, replace it with the implied ranking
  if (any(is.na(rankings))) {
    dn <- dimnames(rankings)
    rankings <- lapply(
      split(rankings, f = seq_len(nrow(rankings))),
      function(x) {
        if (sum(is.na(x)) == 1) x[is.na(x)] <- setdiff(seq_along(x), x)
        return(x)
      }
    )
    rankings <- do.call(rbind, rankings)
    dimnames(rankings) <- dn
  }

  if (!is.null(rho_init)) {
    if (!validate_permutation(rho_init)) stop("rho_init must be a proper permutation")
    if (!(sum(is.na(rho_init)) == 0)) stop("rho_init cannot have missing values")
    if (length(rho_init) != n_items) stop("rho_init must have the same number of items as implied by rankings or preferences")
    rho_init <- matrix(rho_init, ncol = 1)
  }

  # Generate the constraint set
  if (!is.null(preferences) && is.null(constraints)) {
    constraints <- generate_constraints(preferences, n_items)
  } else if (is.null(constraints)) {
    constraints <- list()
  }

  if (is.null(obs_freq)) obs_freq <- rep(1, nrow(rankings))

  logz_list <- prepare_partition_function(logz_estimate, metric, n_items)

  fits <- list(
            run_mcmc_dpmixture(
              rankings = t(rankings),
              obs_freq = obs_freq,
              nmc = nmc,
              constraints = constraints,
              cardinalities = logz_list$cardinalities,
              logz_estimate = logz_list$logz_estimate,
              rho_init = rho_init,
              metric = metric,
              error_model = ifelse(is.null(error_model), "none", error_model),
              Lswap = swap_leap,
              lambda = lambda,
              alpha_max = alpha_max,
              psi = psi,
              psi_init = psi_init,
              leap_size = leap_size,
              alpha_prop_sd = alpha_prop_sd,
              alpha_init = alpha_init,
              alpha_jump = alpha_jump,
              rho_thinning = rho_thinning,
              aug_thinning = aug_thinning,
              clus_thin = clus_thin,
              save_aug = save_aug,
              verbose = verbose,
              kappa_1 = 1.0,
              kappa_2 = 1.0
              )
  )


  if (verbose) {
    print("Metropolis-Hastings algorithm completed. Post-processing data.")
  }

  fit <- tidy_mcmc_dpmixture(
    fits, rho_thinning, rankings, alpha_jump,
    fits[[1]]$max_cluster_index, nmc, aug_thinning,
    n_items, clus_thin)

  fit$save_aug <- save_aug
  fit$rho_thinning <- rho_thinning
  fit$alpha_jump <- alpha_jump
  fit$clus_thin <- clus_thin

  # Add class attribute
  class(fit) <- "BayesMallowsDPMixture"

  return(fit)
}
