# gmdhpmm: GMDH with PMM coefficient estimation

Public code supplement for the Paper 1A methodology manuscript:

> Cumulant-Adaptive Coefficient Estimation in Self-Organized Polynomial Networks

This repository contains the standalone R package `gmdhpmm`, experiment scripts,
and archived result tables for the cumulant-adaptive GMDH-PMM method. Large
licensed raw datasets and monorepo-only manuscript notes are intentionally not
redistributed here.

## Scope

`gmdhpmm` implements MIA-GMDH whose inner coefficient step dispatches between
least squares (LSE), PMM2, and PMM3 using bootstrap-stabilized residual cumulants
of each partial model. It also includes forced baseline modes for ablations:
`ridge-LSE`, `Huber`, and `L1`.

The PMM2/PMM3 Newton-Raphson solvers, moments, and g-factors are delegated to
the school's `EstemPMM` package. The code in this repository provides the
bootstrap cumulant stabilization on small GMDH partitions and the
per-partial-model dispatch inside the self-organized tournament.

## Layout

```text
.
|-- DESCRIPTION, NAMESPACE
|-- R/
|   |-- kg2.R
|   |-- cumulants.R
|   |-- dispatch.R
|   |-- inner_estimate.R
|   |-- external_criterion.R
|   |-- control.R
|   `-- gmdh.R
|-- tests/testthat/
|-- experiments/
|   |-- run_sanity.R
|   |-- run_synthetic.R
|   |-- run_cascade.R
|   |-- run_coverage.R
|   |-- run_realworld.R
|   `-- revision/
`-- datasets/
    `-- external/
```

## Requirements

- R >= 4.1.
- `EstemPMM` for `lm_pmm2` and `lm_pmm3`.
- `testthat`, `pkgload`, and `data.table` for development and tests.
- Optional packages for selected real-data experiments: `MASS`, `survival`,
  `ISLR2`, `boot`, and dataset-specific packages named in the scripts.

## Quick Start

Run from the repository root.

```bash
Rscript -e 'pkgload::load_all(".")'
Rscript -e 'testthat::test_local(".")'
Rscript experiments/run_sanity.R
```

Minimal usage:

```r
pkgload::load_all(".")

X <- matrix(rnorm(600 * 4), 600, 4)
y <- 1 + 1.5 * X[, 1] - X[, 2] + 0.8 * X[, 1] * X[, 2] +
  0.5 * X[, 3]^2 + (rexp(600) - 1) * 2

fit <- gmdh_pmm(X, y, gmdh_pmm_control(B = 200, L_max = 4, F = 6, seed = 7))
fit
yhat <- predict(fit, X)
```

## Verification Commands

The lightweight validation gate is:

```bash
Rscript -e 'testthat::test_local(".")'
Rscript experiments/run_sanity.R
R CMD build .
R CMD check --no-manual --no-build-vignettes gmdhpmm_0.1.0.tar.gz
```

`R CMD check` currently reports one documentation warning because exported
research functions do not yet have `.Rd` help pages. The warning does not affect
the tests or experiment scripts.

The sanity experiment demonstrates two basic behaviours:

- Under Gaussian noise, dispatch falls back overwhelmingly to LSE.
- Under strong right-skew noise, PMM2 is selected for most partial models and
  improves the seeded out-of-sample MSE relative to forced LSE.

## Experiment Scripts

Synthetic and revision experiments are standalone:

```bash
Rscript experiments/run_synthetic.R [R n B]
Rscript experiments/run_cascade.R [R n B]
Rscript experiments/run_coverage.R [R n Bci Bd]
Rscript experiments/revision/e2_threshold_sensitivity.R
Rscript experiments/revision/e3_ablation.R
Rscript experiments/revision/e4_leverage_bootstrap.R
Rscript experiments/revision/e5_inhull_var.R
Rscript experiments/revision/e6_gmm.R
Rscript experiments/revision/e7_significance.R
Rscript experiments/revision/e8_reservetw4_sensitivity.R
```

Real-data screening can be run on built-in/package datasets:

```bash
Rscript experiments/run_realworld.R airquality_ozone 1 10
Rscript experiments/run_realworld.R all 3 20
Rscript experiments/run_realworld.R concrete 3 20
```

Only the UCI Concrete Compressive Strength CSV is shipped in this repository.
Other external datasets listed in `datasets/external/external_candidates.csv`
are documented for provenance, but their raw/prepared files are not
redistributed here. Archived result summaries for those runs are retained under
`experiments/results/`.

Many scripts write CSV/PNG outputs to `experiments/results/`. If you are
running exploratory smoke checks, inspect `git status` before committing so that
new local outputs are not mixed with the archived paper results.

## Archived Results

The repository includes CSV and PNG outputs used during manuscript development.
The main reading is deliberately conservative:

- H2 passes: empirical estimator efficiency tracks the theoretical PMM
  g-factor in PMM-favourable regimes.
- H4 as originally phrased is not supported: OLS does not show exponential bias
  accumulation in the clean cascade experiment.
- The defensible contribution is variance stability, efficient coefficient
  estimation, interval behaviour, and avoiding fragile LSE regimes rather than
  universal prediction-MSE superiority over all robust baselines.

## License and Citation

The package is distributed under GPL-3. See `LICENSE`.

For manuscript-side citation metadata, see `CITATION.cff`.
