# Simulation code for LLnGM Gibbs sampler experiments

This repository contains the simulation code used for the numerical experiments
in the paper:

> Geometric ergodicity of Gibbs samplers for linear latent models with GIG
> variance mixtures

The code reproduces the two simulation studies reported in the manuscript:

- **S1:** mixing comparison across representative GIG parameter regimes.
- **S2:** null-smallness scan obtained by varying the drift parameter `mu`.

## Requirements

The scripts are written in R. They were checked with:

- R 4.4.3
- Matrix 1.7.2
- ngme2 0.9.7
- future 1.58.0
- future.apply 1.20.0
- progressr 0.18.0
- ggplot2 4.0.1
- dplyr 1.1.4
- tidyr 1.3.1

Install the required packages with:

```r
install.packages(c(
  "Matrix",
  "future",
  "future.apply",
  "progressr",
  "ggplot2",
  "dplyr",
  "tidyr"
))
```

The sampler uses `ngme2::rgig()` and `ngme2::ar1()`. Install `ngme2` from the
source used by your project environment if it is not available from your normal
R package repositories.

## Repository Layout

```text
gibbs.R                              Core non-centered Gibbs sampler.
paths.R                              Path helper for standalone/repo-root runs.
test.R                               Small sanity check for the sampler.
experiment.R                         Runs Simulation S1.
make_b1_table.R                      Prints the LaTeX table for S1.
experiment-2.R                       Runs Simulation S2.
Plot-2.R                             Generates the S2 IACT figure.
experiment_summary.csv               Saved S1 summary used in the manuscript.
experiment2_summary_scanmu_4stats.csv Saved S2 summary used for Figure 1.
B2_IACT.png                          Figure 1 in the manuscript.
```

The scripts can be run either from this repository root or from the parent
manuscript repository. The full experiments are computationally heavier than the
sanity check because they use 50,000 MCMC iterations per chain.

## Quick Check

Run a short sampler sanity check:

```sh
Rscript test.R
```

## Reproducing Simulation S1

Run the S1 experiment:

```sh
Rscript experiment.R
```

This writes:

```text
experiment_summary.csv
experiment_results.rds
```

The `.rds` file contains raw chain-level output and is ignored by git because it
can be regenerated. The CSV file is the compact summary used in the manuscript.

To print the LaTeX table from the saved S1 CSV:

```sh
Rscript make_b1_table.R
```

## Reproducing Simulation S2

Run the S2 null-smallness scan:

```sh
Rscript experiment-2.R
```

This writes:

```text
experiment2_summary_scanmu_4stats.csv
experiment2_results_scanmu_4stats.rds
```

The `.rds` file contains raw chain-level output and is ignored by git because it
can be regenerated. The CSV file is used to produce the manuscript figure.

Generate the figure:

```sh
Rscript Plot-2.R
```

This writes:

```text
B2_IACT.png
```

## Reproducibility Notes

The main experiment scripts set random seeds internally:

- S1 uses `set.seed(1)` and chain-specific seeds based on `seed_base`.
- S2 uses `set.seed(2)` and deterministic seeds indexed by scan point and chain.

The GIG sampler is `ngme2::rgig()`. Because `ngme2::rgig()` uses its own
`seed` argument, `gibbs.R` generates deterministic integer seeds from R's seeded
RNG stream and passes them explicitly to `ngme2::rgig(seed = ...)`.

Both experiments use four overdispersed initial states, including fixed vectors
and a random draw from `GIG(1, 1, 1)`.

## Citation

If you use this code, please cite the accompanying paper.
