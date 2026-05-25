# MovieLens Experiment

This folder contains the code to reproduce the MovieLens experiment presented in the paper.
The analysis applies the Dirichlet Process Mixture of Mallows models (DPM3) and a finite
mixture model (FM) to pairwise preference data derived from the MovieLens 1M dataset, and
compares their ability to cluster users and predict held-out preferences.

## Data

The raw data must be downloaded from the GroupLens website:

**https://grouplens.org/datasets/movielens/1m/**

After downloading and extracting the archive, place the following two files in the `data/` folder:

- `ratings.dat` — user ratings
- `users.dat` — user demographic information

The trimmed dataset used in the paper is available directly as `data/data_movielens.RData`,
so the preprocessing step can be skipped if you wish to reproduce the exact results.

## Scripts

The scripts should be run in order:

### `01_preprocess_data.R`

Preprocesses the raw MovieLens data into the pairwise preference format used in the analysis.
The preprocessing pipeline consists of eight steps: selecting the 200 most popular movies,
sampling 300 random assessors, retaining the 50 movies with highest rating entropy, filtering
assessors with insufficient rating diversity, relabelling users and movies, and holding out
one preference per assessor for predictive evaluation. The processed data is saved to
`data/data_movielens.RData`. The paper results use the dataset already available in the
`data/` folder, but this script documents the full preprocessing pipeline.

### `02_run_models.R`

Fits the DPM3 and finite mixture models to the pairwise preference data and computes the
posterior probabilities of correctly predicting the held-out preferences. This script is
computationally intensive and is intended to be run on a server. Results are saved to
`data/result_dpm.RData`, `data/result_fm.RData`, and `data/probabilities.RData`. The burnin
for the co-clustering matrix computation is set conservatively to 100000 and can be adjusted
based on visual inspection of the convergence diagnostics produced by script 04.

### `03_clustering.R`

Applies hierarchical clustering and PAM to the co-clustering matrix obtained from DPM3,
in order to derive a hard partition of the assessors. Produces diagnostic plots including
the dendrogram, PAM elbow, and silhouette plots for both methods. The final partition and
clustering objects are saved to `data/clustering_results.RData`.

### `04_plots.R`

Produces all paper and supplementary figures, including the elbow plot for the finite mixture
model, the DPM3 convergence trace plots, the co-clustering matrix, the posterior probability
boxplots, and the diagonal scatter plot comparing the two methods. Assumes all previous
scripts have been run. Figures are saved to the `plots/` folder.

## Installation

To install the package:

```r
# install.packages("remotes")
remotes::install_github("lorenzo-zuccato/BayesMallows-dpmixture", ref = "dpmixture")
```
