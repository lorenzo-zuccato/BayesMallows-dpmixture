
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- Fork Notice -->

**⚠️ Research Fork of BayesMallows**

This repository extends the official **BayesMallows** package by
implementing a **Dirichlet-process mixture model** for the Mallows
model.  
The purpose of this fork is to allow **automatic discovery of the number
of clusters** in a population of assessors, using a Bayesian
nonparametric approach.

This version contains methods that are **currently undergoing peer
review**. For the stable version of the package, please refer to the
original repository.

# Dirichlet-Process Mixture for the Bayesian Mallows Model

This fork introduces a **Dirichlet-process mixture model** for the
Mallows model, allowing for:

- Automatic detection of the number of clusters among assessors
- Bayesian nonparametric inference
- Integration within the existing BayesMallows framework

## Installation

To install this version of the package use the `remotes` package:

``` r
# install.packages("remotes")
remotes::install_github("lorenzo-zuccato/BayesMallows-dpmixture", ref = "dpmixture")
```

Load the package with

``` r
library(BayesMallowsDPMixture)
```

## Basic Usage Example

Create synthetic data to test the DPMixture:

``` r
set.seed(123)
d1 <- sample_mallows(1:10, 5, 20)
d2 <- sample_mallows(10:1, 2, 20)
data <- rbind(d1,d2)
```

Run the MCMC with `compute_mallows_dpmixture`:

``` r
fit <- compute_mallows_dpmixture(data)
```

Trace plots to assess the convergence of the chain can be generated
using the function `assess_convergence_dpmixture`, which, by default,
displays the trend of the number of non-empty clusters.

``` r
assess_convergence_dpmixture(fit)
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

The same function permits to plot either $\alpha$ or the estimated
cluster proportions $\hat{\tau}$. In these cases, the user must specify
the number of clusters to visualize in the `n` argument. The software
will automatically show the clusters corresponding to the $n$ most
persistent labels in the chain, where the persistence is defined as the
number of samples in which the label is present among the cluster
assignments.

``` r
assess_convergence_dpmixture(fit, parameter = "alpha", n = 5)
```

![](man/figures/README-unnamed-chunk-7-1.png)<!-- -->

``` r
assess_convergence_dpmixture(fit, parameter = "empirical_cluster_probs", n = 5)
```

![](man/figures/README-unnamed-chunk-8-1.png)<!-- -->

Fix the burn-in using the command:

``` r
fit$burnin <- 200
```

Compute the co-clustering matrix with the function
`compute_co_clustering`. This will discard the burn-in and compute the
posterior probability for every two items to belong to the same cluster.

``` r
fit$co_clustering <- compute_co_clustering(fit)
```

Estimate a partition of the assessors after computing the co-clustering
matrix with the function `partition_estimate`. Notice that this method
calls the function `minVI` from the package `mcclust.ext` (see Wade and
Ghahramani ([2015](#ref-wade2015))).

``` r
remotes::install_github("sarawade/mcclust.ext")
```

``` r
library(mcclust.ext)
#> Loading required package: mcclust
#> Loading required package: lpSolve
fit$partition <- partition_estimate(fit)
```

You can visualize the estimated number of non-empty clusters and their
cardinality through the commands:

``` r
max(fit$partition$cl)
#> [1] 2
table(fit$partition$cl)
#> 
#>  1  2 
#> 20 20
```

After a partition has been estimated, you can visualize the
corresponding posterior distributions with the function `plot`, compute
credible intervals with `compute_posterior_intervals` and estimate a
consensus ranking with `compute_consensus`, in the same way as the
original package `BayesMallows`.

``` r
plot(fit, parameter = "alpha")
```

![](man/figures/README-unnamed-chunk-14-1.png)<!-- -->

``` r
compute_posterior_intervals(fit, parameter = "alpha")
#>     cluster parameter  mean median conf_level          hpdi central_interval
#> 1 Cluster 1     alpha 5.359  5.363       95 % [4.581,6.255]    [4.520,6.210]
#> 2 Cluster 2     alpha 1.869  1.847       95 % [1.139,2.584]    [1.153,2.653]
```

# BayesMallows

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BayesMallows)](https://cran.r-project.org/package=BayesMallows)
[![R-CMD-check](https://github.com/ocbe-uio/BayesMallows/workflows/R-CMD-check/badge.svg)](https://github.com/ocbe-uio/BayesMallows/actions)
[![Codecov test
coverage](https://codecov.io/gh/ocbe-uio/BayesMallows/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ocbe-uio/BayesMallows?branch=master)
[![CodeFactor](https://www.codefactor.io/repository/github/ocbe-uio/bayesmallows/badge/develop)](https://www.codefactor.io/repository/github/ocbe-uio/bayesmallows/overview/develop)

This package provides a general framework for analyzing rank and
preference data based on the Bayesian Mallows model first described in
Vitelli et al. ([2018](#ref-vitelli2018)).

## Installation

To install the current release, use

``` r
install.packages("BayesMallows")
```

To install the current development version, use

``` r
# install.packages("remotes")
remotes::install_github("ocbe-uio/BayesMallows")
```

## Basic Usage Example

To get started, load the package with

``` r
library(BayesMallows)
```

The package comes with several example datasets. The simplest one
contains 12 persons’ assessments of the weights of 20 potatoes, either
by visual inspection (`potato_visual`) or by lifting the potatoes and
comparing their relative weights by hand (`potato_weighing`). To fit a
Bayesian Mallows model on the `potato_visual` dataset, we do

``` r
fit <- compute_mallows(potato_visual)
```

Next, we can see a diagnostic plot for the Metropolis-Hastings algorithm
with `assess_convergence()`. The plot below is for the scale parameter,
which measures the variation between the individual rankings.

``` r
assess_convergence(fit)
```

Setting the burnin to 500, we obtain a plot of the posterior
distribution of the scale parameter with:

``` r
plot(fit, burnin = 500)
```

For more examples, please our [R Journal
paper](https://journal.r-project.org/archive/2020/RJ-2020-026/index.html),
and the function documentation. The use of parallel chains are described
in [this
vignette](https://ocbe-uio.github.io/BayesMallows/articles/parallel_chains.html).

## The Bayesian Mallows Model

### Methodology

The BayesMallows package currently implements the complete model
described in Vitelli et al. ([2018](#ref-vitelli2018)), which includes a
large number of distance metrics, handling of missing ranks and pairwise
comparisons, and clustering of users with similar preferences. The
extension to non-transitive pairwise comparisons by Crispino et al.
([2019](#ref-crispino2019)) is also implemented. In addition, the
partition function of the Mallows model can be estimated using the
importance sampling algorithm of Vitelli et al.
([2018](#ref-vitelli2018)) and the asymptotic approximation of Mukherjee
([2016](#ref-mukherjee2016)). For a review of ranking models in general,
see Liu, Crispino, et al. ([2019](#ref-liu2019)). Crispino and
Antoniano-Villalobos ([2022](#ref-crispino2022)) outlines how
informative priors can be used within the model.

Updating of the posterior distribution based on new data, using
sequential Monte Carlo methods, is implemented and described in [a
separate
vignette](https://ocbe-uio.github.io/BayesMallows/articles/SMC-Mallows.html).
The computational algorithms are described in further detail in Stein
([2023](#ref-steinSequentialInferenceMallows2023)).

### Applications

Among the current applications, Liu, Reiner, et al.
([2019](#ref-liu2019b)) applied the Bayesian Mallows model for providing
personalized recommendations based on clicking data, and Barrett and
Crispino ([2018](#ref-barrett2018)) used the model of Crispino et al.
([2019](#ref-crispino2019)) to analyze listeners’ understanding of
music. Eliseussen, Fleischer, and Vitelli
([2022](#ref-eliseussenRankbasedBayesianVariable2022)) presented an
extended model for variable selection in genome-wide transcriptomic
analyses.

### Future Extensions

Plans for future extensions of the package include implementation of a
variational Bayes algorithm for approximation the posterior
distribution. The sequential Monte Carlo algorithms will also be
extended to cover a larger part of the model framework, and we will add
more options for specifications of prior distributions.

## Citation

If using the BayesMallows package in academic work, please cite Sørensen
et al. ([2020](#ref-sorensen2020)), in addition to the relevant
methodological papers.

``` r
citation("BayesMallows")
#> To cite package 'BayesMallows' in publications use:
#> 
#>   Sørensen Ø, Crispino M, Liu Q, Vitelli V (2020). "BayesMallows: An R
#>   Package for the Bayesian Mallows Model." _The R Journal_, *12*(1),
#>   324-342. doi:10.32614/RJ-2020-026
#>   <https://doi.org/10.32614/RJ-2020-026>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     author = {{\O}ystein S{\o}rensen and Marta Crispino and Qinghua Liu and Valeria Vitelli},
#>     doi = {10.32614/RJ-2020-026},
#>     title = {BayesMallows: An R Package for the Bayesian Mallows Model},
#>     journal = {The R Journal},
#>     number = {1},
#>     pages = {324--342},
#>     volume = {12},
#>     year = {2020},
#>   }
```

## Contribution

This is an open source project, and all contributions are welcome. Feel
free to open an
[Issue](https://github.com/ocbe-uio/BayesMallows/issues), a [Pull
Request](https://github.com/ocbe-uio/BayesMallows/pulls), or to e-mail
us.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-barrett2018" class="csl-entry">

Barrett, N., and Marta Crispino. 2018. “The Impact of 3-d Sound
Spatialisation on Listeners’ Understanding of Human Agency in Acousmatic
Music.” *Journal of New Music Research* 47 (5): 399–415.
<https://doi.org/10.1080/09298215.2018.1437187>.

</div>

<div id="ref-crispino2022" class="csl-entry">

Crispino, Marta, and Isadora Antoniano-Villalobos. 2022. “Informative
Priors for the Consensus Ranking in the Bayesian Mallows Model.”
*Bayesian Analysis*, January, 1–24. <https://doi.org/10.1214/22-BA1307>.

</div>

<div id="ref-crispino2019" class="csl-entry">

Crispino, Marta, E. Arjas, V. Vitelli, N. Barrett, and A. Frigessi.
2019. “A Bayesian Mallows Approach to Nontransitive Pair Comparison
Data: How Human Are Sounds?” *The Annals of Applied Statistics* 13 (1):
492–519. <https://doi.org/10.1214/18-aoas1203>.

</div>

<div id="ref-eliseussenRankbasedBayesianVariable2022" class="csl-entry">

Eliseussen, Emilie, Thomas Fleischer, and Valeria Vitelli. 2022.
“Rank-Based Bayesian Variable Selection for Genome-Wide Transcriptomic
Analyses.” *Statistics in Medicine* 41 (23): 4532–53.
<https://doi.org/10.1002/sim.9524>.

</div>

<div id="ref-liu2019" class="csl-entry">

Liu, Q., Marta Crispino, I. Scheel, V. Vitelli, and A. Frigessi. 2019.
“Model-Based Learning from Preference Data.” *Annual Review of
Statistics and Its Application* 6 (1).
<https://doi.org/10.1146/annurev-statistics-031017-100213>.

</div>

<div id="ref-liu2019b" class="csl-entry">

Liu, Q., A. H. Reiner, A. Frigessi, and I. Scheel. 2019. “Diverse
Personalized Recommendations with Uncertainty from Implicit Preference
Data with the Bayesian Mallows Model.” *Knowledge-Based Systems* 186
(December): 104960. <https://doi.org/10.1016/j.knosys.2019.104960>.

</div>

<div id="ref-mukherjee2016" class="csl-entry">

Mukherjee, S. 2016. “Estimation in Exponential Families on
Permutations.” *The Annals of Statistics* 44 (2): 853–75.
<https://doi.org/10.1214/15-aos1389>.

</div>

<div id="ref-sorensen2020" class="csl-entry">

Sørensen, Øystein, Marta Crispino, Qinghua Liu, and Valeria Vitelli.
2020. “BayesMallows: An R Package for the Bayesian Mallows Model.” *The
R Journal* 12 (1): 324–42. <https://doi.org/10.32614/RJ-2020-026>.

</div>

<div id="ref-steinSequentialInferenceMallows2023" class="csl-entry">

Stein, Anja. 2023. “Sequential Inference with the Mallows Model.” PhD
thesis, Lancaster University.

</div>

<div id="ref-vitelli2018" class="csl-entry">

Vitelli, V., Ø. Sørensen, M. Crispino, E. Arjas, and A. Frigessi. 2018.
“Probabilistic Preference Learning with the Mallows Rank Model.”
*Journal of Machine Learning Research* 18 (1): 1–49.
<https://jmlr.org/papers/v18/15-481.html>.

</div>

<div id="ref-wade2015" class="csl-entry">

Wade, Sara, and Zoubin Ghahramani. 2015. “Bayesian Cluster Analysis:
Point Estimation and Credible Balls.” *Bayesian Analysis* 13 (May).
<https://doi.org/10.1214/17-BA1073>.

</div>

</div>
