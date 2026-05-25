# 02_run_models.R
#
# Runs the DPM3 and finite mixture models on the MovieLens dataset and computes
# posterior probabilities of correctly predicting held-out preferences.
# Assumes data/data_movielens.RData has been produced by 01_preprocess_data.R.
# This script is computationally intensive and is intended to be run on a server.

library(BayesMallowsDPMixture)
library(dplyr)
library(mcclust.ext)

# ── Load data ─────────────────────────────────────────────────────────────────

load("data/data_movielens.RData")

# ── Helper function ───────────────────────────────────────────────────────────

compute_probabilities <- function(result, deleted_preferences, ratings_pref, burnin) {
  augmented_data <- result$augmented_data[result$augmented_data$iteration > burnin, ]
  n_iterations   <- length(unique(augmented_data$iteration))
  
  deleted_pref              <- deleted_preferences
  deleted_pref$bottom_item  <- paste0("Item ", deleted_pref$bottom_item)
  deleted_pref$top_item     <- paste0("Item ", deleted_pref$top_item)
  
  probabilities <- inner_join(augmented_data, deleted_pref, by = "assessor") %>%
    filter(item == bottom_item | item == top_item) %>%
    group_by(assessor, iteration) %>%
    summarise(
      count_meeting_criterion = value[item == bottom_item] > value[item == top_item],
      .groups = "drop"
    ) %>%
    group_by(assessor) %>%
    summarise(
      prob = sum(count_meeting_criterion) / n_iterations,
      .groups = "drop"
    )
  
  probabilities$number_preferences <- as.numeric(table(ratings_pref$assessor))
  return(probabilities)
}

# ── DPM3 ──────────────────────────────────────────────────────────────────────

result_dpm <- compute_mallows_dpmixture(
  preferences   = ratings_pref,
  nmc           = 5000,
  clus_thin     = 10,
  rho_thinning  = 10,
  alpha_jump    = 10,
  save_aug      = TRUE,
  aug_thinning  = 10,
  verbose       = TRUE,
  psi           = 0.007,
  alpha_prop_sd = 1,
  leap_size     = 25
)

# Burnin set to 100000; adjust based on trace plot inspection
result_dpm$burnin        <- 1000
result_dpm$co_clustering <- compute_co_clustering(result_dpm)
result_dpm$partition     <- partition_estimate(result_dpm)

prob_dpm <- compute_probabilities(result_dpm, deleted_preferences, ratings_pref, burnin = 1000)

save(result_dpm, file = "data/result_dpm.RData")
cat("DPM3 done\n")

# ── Finite mixture ────────────────────────────────────────────────────────────

result_fm <- compute_mallows_mixtures(
  n_clusters    = 1:15,
  preferences   = ratings_pref,
  nmc           = 5000,
  clus_thin     = 10,
  rho_thinning  = 10,
  alpha_jump    = 10,
  save_aug      = TRUE,
  aug_thinning  = 10,
  verbose       = TRUE,
  include_wcd   = TRUE,
  alpha_prop_sd = 0.5
)

prob_fm <- compute_probabilities(result_fm[[2]], deleted_preferences, ratings_pref, burnin = 1000)

save(result_fm, file = "data/result_fm.RData")
cat("FM done\n")

# ── Save probabilities ────────────────────────────────────────────────────────

save(prob_dpm, prob_fm, file = "data/probabilities.RData")
cat("Probabilities saved\n")