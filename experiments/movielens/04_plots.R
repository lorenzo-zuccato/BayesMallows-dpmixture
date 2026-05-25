# 04_plots.R
#
# Produces all paper figures for the MovieLens experiment.
# Assumes data/result_dpm.RData, data/result_fm.RData, data/probabilities.RData
# and data/clustering_results.RData have been produced by 02 and 03.
library(BayesMallowsDPMixture)
library(ggplot2)
library(dplyr)
options(bitmapType = "cairo")

# ── Load data ─────────────────────────────────────────────────────────────────

load("data/result_dpm.RData")
load("data/result_fm.RData")
load("data/probabilities.RData")

plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

# ── Elbow plot (finite mixture) ───────────────────────────────────────────────

elbow <- plot_elbow(result_fm[2:15], burnin = 100000) +
  theme_bw(base_size = 9)

ggsave(
  filename = file.path(plots_dir, "elbow.png"),
  plot     = elbow,
  width    = 4.7,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── Trace plots (DPM3) ────────────────────────────────────────────────────────

scale_x_iterations <- scale_x_continuous(
  labels = function(x) paste0(x / 1000, "k")
)

c1 <- assess_convergence_dpmixture(result_dpm) +
  scale_x_iterations

ggsave(
  filename = file.path(plots_dir, "trace_n_clusters.png"),
  plot     = c1,
  width    = 2.5,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

c2 <- assess_convergence_dpmixture(result_dpm, parameter = "alpha", n = 10) +
  theme(legend.position = "none") +
  scale_x_iterations

ggsave(
  filename = file.path(plots_dir, "trace_alpha.png"),
  plot     = c2,
  width    = 2.5,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

c3 <- assess_convergence_dpmixture(result_dpm, parameter = "empirical_cluster_probs", n = 10) +
  theme(legend.position = "none") +
  scale_x_iterations

ggsave(
  filename = file.path(plots_dir, "trace_cluster_probs.png"),
  plot     = c3,
  width    = 2.5,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── Co-clustering matrix (DPM3) ───────────────────────────────────────────────

co_clust <- plot(result_dpm, parameter = "co_clustering") +
  theme_bw(base_size = 9) +
  coord_fixed() +
  theme(
    legend.position   = "bottom",
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.text       = element_text(size = 6),
    legend.title      = element_blank(),
    axis.text         = element_text(size = 7),
    axis.ticks        = element_line(linewidth = 0.3),
    legend.margin     = margin(t = -7.5)
  ) +
  scale_fill_viridis_c(option = "viridis")

ggsave(
  filename = file.path(plots_dir, "co_clust.png"),
  plot     = co_clust,
  width    = 2.5,
  height   = 2.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── Probability boxplots ──────────────────────────────────────────────────────

custom_breaks <- c(0, 30, 60, 100, 150, 200, 300, max(prob_dpm$number_preferences))
custom_labels <- sprintf("[%d - %d]",
                         custom_breaks[1:(length(custom_breaks) - 1)],
                         custom_breaks[2:length(custom_breaks)] - 1)

make_boxplot <- function(df) {
  df$preferences_range <- cut(df$number_preferences, breaks = custom_breaks,
                              labels = as.character(seq(1, length(custom_breaks) - 1)))
  df_long <- df

  boxplot_plot <- ggplot(df_long, aes(x = factor(preferences_range), y = prob)) +
    geom_boxplot(fill = "lightgrey") +
    labs(x = "Number of preferences",
         y = "Posterior probability\nof correct preference prediction") +
    scale_x_discrete(labels = custom_labels) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
    )

  histogram_plot <- ggplot(df_long, aes(x = prob)) +
    geom_histogram(bins = 20, fill = "darkgrey", color = "black") +
    theme_bw(base_size = 9) +
    theme(
      axis.text.y  = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    coord_flip()

  patchwork::wrap_plots(boxplot_plot, histogram_plot, ncol = 2, widths = c(8, 2))
}

ggsave(
  filename = file.path(plots_dir, "prob_fm.png"),
  plot     = make_boxplot(prob_fm),
  width    = 3.5,
  height   = 3.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

ggsave(
  filename = file.path(plots_dir, "prob_dpm.png"),
  plot     = make_boxplot(prob_dpm),
  width    = 3.5,
  height   = 3.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)

# ── Diagonal scatter plot ─────────────────────────────────────────────────────

your_data <- data.frame(
  prob_1      = prob_fm$prob,
  prob_2      = prob_dpm$prob,
  n_pref_max  = prob_dpm$number_preferences
)

diagonal_plot <- ggplot(your_data, aes(x = prob_1, y = prob_2, color = n_pref_max)) +
  geom_point(size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_viridis_c(
    name   = "No. preferences",
    option = "viridis",
    trans  = "sqrt",
    breaks = c(1, 25, 80, 180, 300),
    labels = c("1", "25", "80", "180", "300+")
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position   = "bottom",
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7),
    legend.margin     = margin(t = -5),
    aspect.ratio      = 1
  ) +
  labs(
    x = "Posterior probability — finite mixture",
    y = "Posterior probability — DPM3"
  )

ggsave(
  filename = file.path(plots_dir, "diagonal_plot.png"),
  plot     = diagonal_plot,
  width    = 3.5,
  height   = 3.5,
  units    = "in",
  dpi      = 600,
  device   = "png"
)
