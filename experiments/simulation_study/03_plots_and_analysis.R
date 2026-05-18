library(patchwork)
library(BayesMallowsDPMixture)
library(mclust)
library(ggplot2)
options(bitmapType = "cairo")

# ── Settings ──────────────────────────────────────────────────────────────────

# Change data_type to "top8" or "pref30" as needed
data_type     <- "pref30"
true_clusters <- c(rep(1, 50), rep(2, 25), rep(3, 25))

# ── Output directory ──────────────────────────────────────────────────────────

plots_dir <- file.path("plots", data_type)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ── Accumulators ──────────────────────────────────────────────────────────────

dpm_wins <- 0
fm_wins  <- 0
ties     <- 0

all_ari_dpm  <- numeric(30)
all_ari_fm   <- numeric(30)
all_n_differ <- numeric(30)
all_n_clust  <- numeric(30)

ari_dpm_3  <- c()
ari_fm_3   <- c()
n_differ_3 <- c()

# ── Main loop ─────────────────────────────────────────────────────────────────

for (i in 1:30) {
  load(file.path("results", data_type, paste0("result_", data_type, "_", i, ".RData")))
  
  fm <- result$fm
  dpm <- result$dpm
  rm(result)
  
  # ── Plots ───────────────────────────────────────────────────────────────────
  
  co_clust <- plot(dpm, parameter = "co_clustering") +
    theme_bw(base_size = 11) +
    coord_fixed() +
    theme(
      legend.position   = "bottom",
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(0.6, "cm"),
      legend.text       = element_text(size = 8),
      legend.title      = element_blank(),
      axis.text         = element_text(size = 8),
      axis.ticks        = element_line(linewidth = 0.3),
      legend.margin     = margin(t = -15),
    ) +
    scale_fill_viridis_c(option = "viridis")
  
  elbow <- plot_elbow(fm, burnin = 20000) +
    theme_bw(base_size = 11)
  
  alpha_trace <- assess_convergence_dpmixture(dpm, parameter = "alpha", n = 10) +
    theme_bw(base_size = 11)
  
  cluster_probs_trace <- assess_convergence_dpmixture(
    dpm, parameter = "empirical_cluster_probs", n = 10) +
    theme_bw(base_size = 11)
  
  combined <- (co_clust | elbow) / (alpha_trace | cluster_probs_trace) +
    plot_annotation(title = paste0(data_type, " — dataset ", i))
  
  ggsave(
    filename = file.path(plots_dir, paste0(data_type, "_", i, "_diagnostics.png")),
    plot     = combined,
    width    = 12,
    height   = 10,
    units    = "in",
    dpi      = 150,
    device   = "png"
  )
  
  # ── Partition comparison ────────────────────────────────────────────────────
  
  dpm_clusters   <- dpm$partition$cl
  n_dpm_clusters <- length(unique(dpm_clusters))
  
  fm_partition <- assign_cluster(fm[[3]], burnin = 2000, soft = FALSE, expand = FALSE)
  fm_clusters  <- fm_partition$map_cluster[order(fm_partition$assessor)]
  
  n_differ    <- length(classError(fm_clusters, dpm_clusters)$misclassified)
  ari_dpm     <- adjustedRandIndex(dpm_clusters, true_clusters)
  ari_fm      <- adjustedRandIndex(fm_clusters,  true_clusters)
  n_wrong_dpm <- length(classError(dpm_clusters, true_clusters)$misclassified)
  n_wrong_fm  <- length(classError(fm_clusters,  true_clusters)$misclassified)
  
  closer <- if (ari_dpm > ari_fm) {
    dpm_wins <- dpm_wins + 1; "DPM3"
  } else if (ari_fm > ari_dpm) {
    fm_wins <- fm_wins + 1; "FM"
  } else {
    ties <- ties + 1; "Tie"
  }
  
  all_ari_dpm[i]  <- ari_dpm
  all_ari_fm[i]   <- ari_fm
  all_n_differ[i] <- n_differ
  all_n_clust[i]  <- n_dpm_clusters
  
  if (n_dpm_clusters == 3) {
    ari_dpm_3  <- c(ari_dpm_3,  ari_dpm)
    ari_fm_3   <- c(ari_fm_3,   ari_fm)
    n_differ_3 <- c(n_differ_3, n_differ)
  }
  
  cat("Dataset", i, "\n")
  cat("  DPM3 number of clusters:", n_dpm_clusters, "\n")
  cat("  Assessors in different clusters (DPM3 vs FM):", n_differ,
      "out of", length(dpm_clusters), "\n")
  cat("  Misclassified vs truth — DPM3:", n_wrong_dpm, "| FM:", n_wrong_fm, "\n")
  cat("  ARI vs truth            — DPM3:", round(ari_dpm, 3), "| FM:", round(ari_fm, 3), "\n")
  cat("  Closer to truth:", closer, "\n\n")
}

# ── Summary ───────────────────────────────────────────────────────────────────

cat("=== Summary over all 30 datasets ===\n")
cat("  DPM3 closer to truth:", dpm_wins, "\n")
cat("  FM closer to truth:  ", fm_wins,  "\n")
cat("  Ties:                ", ties,     "\n")
cat("  Avg ARI — DPM3:", round(mean(all_ari_dpm), 3),
    "| FM:", round(mean(all_ari_fm), 3), "\n")
cat("  Avg assessors differing (DPM3 vs FM):", round(mean(all_n_differ), 1), "\n")
cat("  DPM3 found 3 clusters in", sum(all_n_clust == 3), "out of 30 datasets\n")

cat("\n=== Summary restricted to datasets where DPM3 finds 3 clusters (n =",
    length(ari_dpm_3), ") ===\n")
cat("  Avg ARI — DPM3:", round(mean(ari_dpm_3), 3),
    "| FM:", round(mean(ari_fm_3), 3), "\n")
cat("  Avg assessors differing (DPM3 vs FM):", round(mean(n_differ_3), 1), "\n")
cat("  DPM3 closer to truth:", sum(ari_dpm_3 > ari_fm_3), "\n")
cat("  FM closer to truth:  ", sum(ari_fm_3 > ari_dpm_3), "\n")
cat("  Ties:                ", sum(ari_dpm_3 == ari_fm_3), "\n")