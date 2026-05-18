# Simulation Study

This folder contains the code and data to reproduce the simulation study in the paper.
The study compares the performance of the Dirichlet Process Mixture of Mallows models (DPM3)
and finite mixture models (FM) in recovering a known cluster structure from partial ranking
and pairwise preference data.

## Data

Two observation formats are considered:

- **top8**: aach assessor provides a top-k ranking of 30 items, where k averages 8
- **pref30**: each assessor provides on average 30 pairwise preferences among 30 items

In both cases, 100 assessors are simulated from a mixture of three Mallows models with the
following parameters:

| Cluster | Consensus ranking $\rho$       | Precision $\alpha$ | Size |
|---------|--------------------------------|--------------------|------|
| 1       | $(1, \dots, 30)$               | 3                  | 50   |
| 2       | $(30, \dots, 1)$               | 2                  | 25   |
| 3       | $(16, \dots, 30, 1, \dots, 15)$| 5                  | 25   |

30 independent datasets are generated for each format.

The `data/` folder contains:

- `data_top8_paper.RData` — the exact top8 datasets used in the paper
- `data_pref30_paper.RData` — the exact pref30 datasets used in the paper

Each file contains an R list with two elements:
- `$parameters` — the true parameters used in data generation
- `$datasets` — a list of 30 datasets

## Scripts

The scripts should be run in order:

### `01_generate_data.R`

Generates 30 independent datasets for both the top8 and pref30 scenarios and saves them to
`data/data_top8_new.RData` and `data/data_pref30_new.RData`. The paper results use the
datasets in `data/data_top8_paper.RData` and `data/data_pref30_paper.RData`, but this script
shows how they were generated and can be used to produce new datasets of the same type.

### `02_run_models.R`

Runs the DPM3 and FM models for a single dataset. Takes two command line arguments:

- `<index>`: integer from 1 to 30
- `<data_type>`: either `top8` or `pref30`

Example:
```bash
Rscript 02_run_models.R 1 top8
```

Run this for all combinations of index and data type before proceeding to script 03.
To run all 30 datasets for both scenarios sequentially:

```bash
for i in $(seq 1 30); do Rscript 02_run_models.R $i top8; done
for i in $(seq 1 30); do Rscript 02_run_models.R $i pref30; done
```

Results are saved to `results/top8/` and `results/pref30/`. The burnin for the co-clustering
matrix computation is set conservatively to 25000 and can be adjusted at the top of the script
based on visual inspection of the convergence diagnostics produced by script 03.

### `03_plots_and_analysis.R`

Produces diagnostic and results plots for all 30 datasets and prints a summary of the
partition comparison between DPM3 and FM. Set `data_type` at the top of the script to
either `"top8"` or `"pref30"` before running. Assumes all 30 results for the chosen
scenario have been produced by script 02.

For each dataset, a combined diagnostic figure is saved to `plots/<data_type>/` containing:
- Co-clustering matrix
- Elbow plot
- Alpha trace plot
- Empirical cluster probabilities trace plot

At the end of the script, a summary is printed comparing the partition recovery performance
of DPM3 and FM across all 30 datasets, both overall and restricted to datasets where DPM3
correctly identifies 3 clusters.

## Installation

To install the package:

```r
# install.packages("remotes")
remotes::install_github("lorenzo-zuccato/BayesMallows-dpmixture", ref = "dpmixture")
```