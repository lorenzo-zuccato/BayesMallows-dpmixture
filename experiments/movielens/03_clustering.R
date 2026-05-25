# 03_clustering.R
#
# Applies hierarchical clustering and PAM to the co-clustering matrix from
# result_dpm, produces diagnostic plots, and saves the final partition.
# Assumes data/result_dpm.RData has been produced by 02_run_models.R.

library(ggplot2)
library(ggdendro)
library(cluster)
library(mclust)
options(bitmapType = "cairo")

# ── Load data ─────────────────────────────────────────────────────────────────

load("data/result_dpm.RData")

plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

# ── Distance matrix ───────────────────────────────────────────────────────────

similarity_matrix <- result_dpm$co_clustering
distance_matrix   <- as.dist(1 - similarity_matrix)

# ── Hierarchical clustering dendrogram ────────────────────────────────────────

hclust_result <- hclust(distance_matrix, method = "complete")
dendro_data   <- dendro_data(hclust_result)

dendrogram <- ggplot() +
  geom_segment(
    data = dendro_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank()
  ) +
  xlab("Assessors") +
  ylab("Height")

ggsave(
  filename = file.path(plots_dir, "dendrogram.png"),
  plot     = dendrogram,
  width    = 4.7,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── PAM elbow ─────────────────────────────────────────────────────────────────

k_values <- 1:10
withinss <- numeric(length(k_values))
for (k in k_values) {
  withinss[k] <- pam(distance_matrix, k = k)$objective["swap"]
}

elbow_pam <- ggplot(
  data.frame(k = k_values, withinss = withinss),
  aes(x = factor(k), y = withinss, group = 1)
) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 9) +
  xlab("Number of clusters") +
  ylab("PAM objective")

ggsave(
  filename = file.path(plots_dir, "elbow_pam.png"),
  plot     = elbow_pam,
  width    = 2.5,
  height   = 3.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── Silhouette plot (PAM) ─────────────────────────────────────────────────────

k_values_sil <- 2:10
sil_pam      <- numeric(length(k_values_sil))
for (k in k_values_sil) {
  sil_pam[k - 1] <- mean(silhouette(pam(distance_matrix, k = k))[, 3])
}

silhouette_pam <- ggplot(
  data.frame(k = k_values_sil, silhouette = sil_pam),
  aes(x = factor(k), y = silhouette, group = 1)
) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 9) +
  xlab("Number of clusters") +
  ylab("Average silhouette width") +
  ylim(0, NA)

ggsave(
  filename = file.path(plots_dir, "silhouette_pam.png"),
  plot     = silhouette_pam,
  width    = 3.5,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── Silhouette plot (hierarchical) ────────────────────────────────────────────

sil_hclust <- numeric(length(k_values_sil))
for (k in k_values_sil) {
  clusters           <- cutree(hclust_result, k = k)
  sil_hclust[k - 1] <- mean(silhouette(clusters, distance_matrix)[, 3])
}

silhouette_hclust <- ggplot(
  data.frame(k = k_values_sil, silhouette = sil_hclust),
  aes(x = factor(k), y = silhouette, group = 1)
) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 9) +
  xlab("Number of clusters") +
  ylab("Average silhouette width") +
  ylim(0, NA)

ggsave(
  filename = file.path(plots_dir, "silhouette_hclust.png"),
  plot     = silhouette_hclust,
  width    = 3.5,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── Partition comparison ──────────────────────────────────────────────────────

pam_clusters    <- pam(distance_matrix, k = 2)$clustering
hclust_clusters <- cutree(hclust_result, k = 2)

n_differ <- length(classError(pam_clusters, hclust_clusters)$misclassified)
cat("Assessors in different clusters (PAM vs hierarchical):", n_differ, "\n")
print(table(PAM = pam_clusters, Hierarchical = hclust_clusters))

# ── Save partition ────────────────────────────────────────────────────────────

# We use the hierarchical clustering partition in the subsequent analysis
partition <- hclust_clusters
save(partition, hclust_result, pam_clusters, file = "data/clustering_results.RData")