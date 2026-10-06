# gmdhpmm: GMDH with PMM coefficient estimation

Public code supplement for the Paper 1A methodology manuscript:

> Cumulant-Adaptive Coefficient Estimation in Self-Organized Polynomial Networks

This repository contains the standalone R package `gmdhpmm`, experiment scripts,
and archived result tables for the cumulant-adaptive GMDH-PMM method. Large
licensed raw datasets and monorepo-only manuscript notes are intentionally not
redistributed here.

## Revision R1 (branch `r1`, 2026-10-06)

The tag `kybernetika-7414-submitted` marks the code as submitted. This branch carries the
corrections prepared for the revision. The package is `gmdhpmm` 0.2.0, the version on
which every R1 number was computed.

- All estimator arms of the cascade and real-data experiments share one inner
  training/validation partition per replication. Earlier, each arm drew its own, which
  mixed estimator effects with split noise.
- Synthetic noise is standardized with population moments, not sample moments. Sample
  centring made the noise orthogonal to the intercept and halved the measured bimodal ARE.
- Summaries report medians with paired tests. Coverage is analysed at the seed level.

Run every script from the repository root. Pre-fix outputs are kept next to the `_v2` ones.

| Manuscript item | Script | Output |
|---|---|---|
| Tab. synthetic, Fig. 1 | `Rscript experiments/run_synthetic.R 1000 2000 500 _v2` | `experiments/results/synthetic_h2_v2.csv` |
| Tab. cascade | `Rscript experiments/run_cascade.R 200 2000 200` | `experiments/results/cascade_recovery_v2{,_raw}.csv` |
| Tab. ablation | `Rscript experiments/revision/e3_ablation.R 50 2000 200 _v2` | `experiments/revision/e3_ablation_v2.csv` |
| Tab. realdata | `Rscript experiments/run_realworld.R <id> 30 200 _v2` for `islr_credit_balance`, `concrete`, `airquality_ozone`, `mass_cars93_mpg`, `survival_solder_skips`, then `Rscript experiments/revision/e11_realdata_v2.R` | `experiments/results/realworld_c5_<id>_v2_*.csv`, `experiments/revision/e11_realdata_v2.csv` |
| Bootstrap-count check | `Rscript experiments/run_realworld.R islr_credit_balance 30 {80,500} _v2_B{80,500}` | `experiments/results/realworld_c5_islr_credit_balance_v2_B*.csv` |
| Coverage (seed level) | `Rscript experiments/revision/e10_coverage_cluster.R` | `experiments/revision/e10_coverage_cluster.csv` |
| PMM3 ARE vs n (supplement) | `Rscript experiments/revision/e9_pmm3_nsweep.R 400` | `experiments/revision/e9_pmm3_nsweep.csv` |
| Tab. synthmech, VAR audit | `Rscript experiments/run_synthetic_var_break.R` (unchanged; fixed split for all arms) | `experiments/results/synthetic_var_audit_*.csv` |

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
