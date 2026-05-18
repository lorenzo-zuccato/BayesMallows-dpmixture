# Usage: Rscript 02_run_models.R <dataset_index> <data_type>
# Example: Rscript 02_run_models.R 1 top8
# <dataset_index>: integer from 1 to 30
# <data_type>: either 'top8' or 'pref30'

library(BayesMallowsDPMixture)
library(mcclust.ext)

# ── Parse arguments ───────────────────────────────────────────────────────────

args      <- commandArgs(trailingOnly = TRUE)
idx       <- as.integer(args[1])
data_type <- args[2]

stopifnot(idx >= 1 && idx <= 30)
stopifnot(data_type %in% c("top8", "pref30"))

# ── Load data ─────────────────────────────────────────────────────────────────

if (data_type == "top8") {
  load("data/data_top8_paper.RData")
  dataset <- data_top8_paper$datasets[[idx]]
} else {
  load("data/data_pref30_paper.RData")
  dataset <- data_pref30_paper$datasets[[idx]]
}

# ── Run models ────────────────────────────────────────────────────────────────

result <- list()

if (data_type == "top8") {
  result$dpm <- compute_mallows_dpmixture(
    rankings     = dataset,
    nmc          = 50000,
    clus_thin    = 10,
    rho_thinning = 10,
    alpha_jump   = 10,
    save_aug     = FALSE,
    verbose      = TRUE,
    psi          = 0.025
  )
  result$fm <- compute_mallows_mixtures(
    n_clusters   = 2:6,
    rankings     = dataset,
    nmc          = 50000,
    clus_thin    = 10,
    rho_thinning = 10,
    alpha_jump   = 10,
    save_aug     = FALSE,
    verbose      = TRUE,
    include_wcd  = TRUE
  )
} else {
  result$dpm <- compute_mallows_dpmixture(
    preferences  = dataset,
    nmc          = 10000,
    clus_thin    = 10,
    rho_thinning = 10,
    alpha_jump   = 10,
    save_aug     = FALSE,
    verbose      = TRUE,
    psi          = 0.025
  )
  result$fm <- compute_mallows_mixtures(
    n_clusters   = 2:6,
    preferences  = dataset,
    nmc          = 50000,
    clus_thin    = 10,
    rho_thinning = 10,
    alpha_jump   = 10,
    save_aug     = FALSE,
    verbose      = TRUE,
    include_wcd  = TRUE
  )
}
# ── Compute co-clustering matrix and partition ────────────────────────────────

# Burnin set conservatively to 25000; adjust based on trace plot inspection
result$dpm$burnin    <- 25000
result$dpm$co_clustering <- compute_co_clustering(result$dpm)
result$dpm$partition     <- partition_estimate(result$dpm)

# ── Save results ──────────────────────────────────────────────────────────────

results_dir <- file.path("results", data_type)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
save(result, file = file.path(results_dir, paste0("result_", data_type, "_", idx, ".RData")))