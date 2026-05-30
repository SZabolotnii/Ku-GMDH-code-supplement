# gmdhpmm — GMDH with PMM coefficient estimation

Reference implementation for **Paper 1** of the GMDH↔Kunchenko program
(`../algorithm-spec.md`, `../problem-statement.md`). MIA-GMDH whose inner
coefficient step dispatches automatically between **LSE**, **PMM2** and
**PMM3**, keyed on **bootstrap-stabilized residual cumulants** of each partial
model. It also includes forced baseline modes for ablations: `ridge-LSE`,
`Huber`, and `L1`.

The hard PMM estimators (Newton–Raphson `lm_pmm2`/`lm_pmm3`, moments, g-factors)
are delegated to the school's **`EstemPMM`** package. The novel pieces here are
the bootstrap cumulant stabilization on small GMDH partitions (contribution C2)
and the per-partial-model dispatch inside the tournament (C1).

## Layout

```
code/
├── DESCRIPTION, NAMESPACE          # R package metadata (Imports: EstemPMM, stats, ...)
├── R/
│   ├── kg2.R                 # KG-2 partial model (features, predict, coef extraction)
│   ├── cumulants.R           # sample + bootstrap-stabilized cumulant diagnostics (C2)
│   ├── dispatch.R            # LSE / PMM2 / PMM3 selection with significance gate (4.1)
│   ├── inner_estimate.R      # INNER_ESTIMATE: warm-start LSE -> diag -> dispatch -> PMM refit/baselines
│   ├── external_criterion.R  # MSE (default) + experimental PMM-loss (8)
│   ├── control.R             # gmdh_pmm_control() tuning parameters (defaults from 9)
│   └── gmdh.R                # MIA tournament, predict(), print()
├── tests/testthat/           # unit + end-to-end tests (testthat edition 3)
└── experiments/run_sanity.R  # Level-1 sanity experiment
```

## Requirements

- R ≥ 4.1, package **EstemPMM** (provides `lm_pmm2`, `lm_pmm3`).
- `Suggests`: `testthat` (≥ 3.0). Dev workflow uses `pkgload`.

## Usage

```r
pkgload::load_all("paper-1-gmdh-pmm/code")   # or install the package

X <- matrix(rnorm(600 * 4), 600, 4)
y <- 1 + 1.5 * X[,1] - X[,2] + 0.8 * X[,1] * X[,2] + 0.5 * X[,3]^2 + (rexp(600) - 1) * 2
fit <- gmdh_pmm(X, y, gmdh_pmm_control(B = 200, L_max = 4, F = 6, seed = 7))
fit                          # summary: depth, best CR, method counts
yhat <- predict(fit, X)
```

## Tests & sanity

```bash
cd paper-1-gmdh-pmm/code
Rscript -e 'testthat::test_local(".")'         # unit + integration suite
Rscript experiments/run_sanity.R               # Level-1 sanity experiment
Rscript experiments/run_synthetic.R [R n B]    # Level-2 H2 benchmark (default 200 500 200)
Rscript experiments/run_cascade.R [R n B]      # GMDH-level recovery + H4 (default 40 800 150)
Rscript experiments/run_coverage.R [R n Bci Bd] # H1 synthetic coverage pilot (default 60 500 80 80)
Rscript experiments/run_realworld.R [candidate|all|domain|extended|external] [R B] # C5 candidate-stage full GMDH comparison
Rscript experiments/run_fremtpl2_sensitivity.R [R B] # freMTPL2 auto/relaxed/forced PMM2 sensitivity
Rscript experiments/run_fremtpl2_reserve_criterion.R [R B] # freMTPL2 reserve-aware selection
Rscript experiments/run_fremtpl2_tail_sweep.R [R B] # freMTPL2 TailReserve weight sweep + decile calibration
Rscript experiments/run_gasturbine_timeaware.R [CO|NOX] [B] # Gas Turbine chronological split + spike calibration
Rscript experiments/run_gasturbine_arx_spike.R [CO|NOX] [B] [tail_weights] # Gas Turbine ARX lags + reserve/spike-aware sweep
Rscript experiments/run_gasturbine_arx_ablation.R [CO|NOX] [B] # Gas Turbine compact ARX feature ablation
Rscript experiments/run_gasturbine_interval_coverage.R [CO|NOX] [B] # Gas Turbine validation-calibrated interval coverage
```

The sanity experiment demonstrates the two behaviours that anchor Paper 1's
hypotheses:

- **H3 (safe fallback):** under Gaussian noise the dispatch is overwhelmingly
  LSE (≈ all layer-0 models; rare harmless deep-layer PMM2 where intermediate
  residuals drift from Gaussian).
- **Dispatch + utility:** under strong right-skew noise PMM2 is selected for the
  bulk of partial models, and GMDH-PMM beats an LSE-forced tournament out of
  sample (≈ 12 % lower test MSE in the seeded run).

`run_synthetic.R` validates **H2** (efficiency scales with the g-factor) at the
estimator level (R = 200, n = 500): empirical efficiency tracks theory with
Pearson **rho = 0.94**, PMM cells show **39–68 %** variance reduction
(Wilcoxon p <= 1e-19), and the leptokurtic-symmetric decoys (Laplace,
contaminated normal) correctly stay on LSE. The same run now reports Huber/L1
baselines: in PMM-favorable cells PMM beats Huber by **2.88–4.91x** and L1 by
**7.46–46.56x** at coefficient-MSE level. Results land in
`experiments/results/` (`synthetic_h2.csv`, `are_emp_vs_theory.png`). **H2: PASS.**

`run_cascade.R` runs a full tournament to depth 4 on a known nonlinear target
and tracks recovery RMSE / bias / variance by depth under Gaussian and skew
noise. Arms are PMM-auto, forced LSE, forced `ridge-LSE` (`lambda = 1e-8`,
intercept unpenalized), Huber, and L1.
**H4: PARTIAL but cleaner.** Coefficient bias stays ≈ 0 in all arms (OLS is
unbiased — no exponential bias accumulation), so H4 as originally stated is not
observed. Under skew at depth 4, PMM gives the best balance in the default run:
recovery RMSE **0.131** vs Huber **0.137**, L1 **0.208**, LSE **0.271**,
ridge-LSE **0.393**; RMSE SD **0.064** vs **0.069**, **0.165**, **0.397**,
**1.288**. Huber is close in RMSE but carries a systematic bias (`-0.148`),
whereas PMM is near zero (`-0.0017`). Takeaway for the paper: PMM-GMDH's value
is **parameter efficiency (H2), cascade variance stability under skew, and
inferential coverage (H1)**, not generic conditional-mean MSE superiority.

`run_coverage.R` runs the synthetic H1 pilot for coefficient intervals. It
stores `coverage_h1.csv`, `coverage_h1_raw.csv`, and `coverage_h1.png`. Current
reading: the strong "LSE collapses to 70% coverage" claim is not supported in
the clean KG-2 synthetic setup, but PMM bootstrap-percentile intervals keep
near-nominal coefficient-wise coverage with much shorter widths. Examples:
skew PMM2-favorable `exp` has PMM coverage **0.939** with width **0.125** vs
LSE classical **0.928** with width **0.171**; platykurtic `uniform` has PMM
coverage **0.947** with width **0.104** vs LSE classical **0.961** with width
**0.173**. Ukrainian analysis: `../coverage-h1-report.md`.

`run_realworld.R` runs C5 full GMDH comparisons. `all` is retained as the
built-in shortlist (`airquality`, `iris`, `ToothGrowth`); `domain` runs Acme
returns, credit balance, soldering skips, motorcycle acceleration, and Cars93
MPG; `external` reads `shared/datasets/external/external_candidates.csv`.
Current reading: no local candidate is a publication-grade flagship.
The best domain near-tie is `islr_credit_balance`: PMM mean NRMSE **0.279** vs
LSE **0.281**. Built-in `airquality_ozone` screened as PMM2-favorable, but full
GMDH mean NRMSE is PMM **0.665** vs LSE **0.636**. The first external candidate
is `freMTPL2` insurance severity: raw claim amounts screen as PMM2 with
expected ARE **8.105**; candidate-stage GMDH shows PMM improving over LSE/ridge
on NRMSE while Huber/L1 trade lower MAE for strong negative severity bias.
`run_fremtpl2_sensitivity.R` then checks auto vs relaxed vs forced PMM2 on raw
severity. Current reading: forced/relaxed PMM2 does not unlock a large
prediction win; robust baselines still minimize MAE/NRMSE but under-reserve
aggregate and top-decile claims heavily, so the next experiment should use a
reserve-aware criterion. Ukrainian analyses: `../realworld-c5-report.md`,
`../domain-c5-report.md`, `../external-c5-fremtpl2-report.md`, and
`../fremtpl2-sensitivity-report.md`.

`run_fremtpl2_reserve_criterion.R` compares classical MSE node selection with
aggregate and tail reserve-aware selection. Current reading: reserve-aware
criteria substantially improve the PMM/LSE family relative to MSE selection
(e.g. PMM-auto NRMSE **1.252 -> 1.045**, tail reserve error **-0.854 -> -0.828**),
but they do not yet establish a PMM-specific win over LSE. Robust Huber/L1
remain best by MAE/NRMSE while heavily under-reserving the tail.

`run_fremtpl2_tail_sweep.R` sweeps `tail_weight` in the TailReserve criterion
and writes decile calibration tables. Current reading: increasing
`tail_weight` from 1 to 16 improves PMM-auto top-decile reserve error only from
**-0.831** to **-0.821**, while LSE shows a similar response (**-0.830** to
**-0.825**). Even the best PMM-auto setup predicts about **3498** average claim
severity in the top decile against observed **20842**. So freMTPL2 is useful as
a cautionary actuarial case, but not yet a flagship PMM-specific win.

After the second DeepResearch pass, `run_realworld.R` also accepts external
UCI Gas Turbine manifest ids prepared by `shared/datasets/prepare_gasturbine.R`.
Current smoke results (`R=10`, `B=30`): `gas_turbine_co_2015_raw` is the first
industrial-emissions survivor, with PMM NRMSE **0.542** vs LSE **0.549** and
Huber **0.540**, while PMM2 is selected in **84.4%** of internal nodes.
NOx is a weaker control target (LSE **0.485**, PMM **0.488**).

`run_gasturbine_timeaware.R` uses a chronological split for the industrial
emissions survivor: train 2011-2013, internal validation 2014, external test
2015. Current CO result: PMM is much better than LSE/ridge-LSE by NRMSE
(**0.767** vs **1.132**) and selects PMM2 in all selected PMM nodes, but does
not improve spike calibration (top-10 reserve error **-0.289** vs LSE
**-0.221**, L1 **-0.140**). This keeps CO as a survivor, not yet a flagship.

`run_gasturbine_arx_spike.R` follows the same chronological split but adds
ARX-like features (`CO_lag1`, `CO_lag2`, and lag-1 operational inputs) and
sweeps reserve/spike-aware validation criteria. Current CO run (`B=30`,
`tail_weight = 2,4,8`): ARX residuals still dispatch to PMM2, while validation
NRMSE improves from **0.768** (static warm start) to **0.573**. PMM with
ReserveTW2/4/8 improves from static PMM **0.767** NRMSE and top-10 reserve
**-0.289** to **0.563** and **-0.048**. SpikeTW4/8 gives near-zero PMM
top-10 reserve (**0.013**) but worse NRMSE (**0.688**). Strong case study,
not yet a clean flagship because L1/Huber remain competitive.

`run_gasturbine_arx_ablation.R` checks whether the ARX gain needs the full
20-feature lag matrix. Current CO run (`B=30`): under ReserveTW4, compact
`target_lags` (static inputs + `CO_lag1/CO_lag2`, 11 features) gives the best
PMM NRMSE **0.556**, slightly better than `full_arx` PMM **0.563** and much
better than `static` PMM **0.706** / `input_lags` PMM **0.610**. This supports
a compact dynamic soft-sensor story: target memory carries most of the PMM
gain, while full input lags mainly improve tail reserve calibration.

`run_gasturbine_interval_coverage.R` calibrates prediction intervals on 2014
validation residuals and tests them on the external 2015 holdout. Current CO
run (`B=30`, target_lags): with ReserveTW4 and absolute-residual 95% intervals,
PMM reaches near-nominal test coverage **0.953** with mean width **3.552**,
versus LSE **0.963** with width **4.168**. PMM is more efficient overall, but
LSE covers the extreme top-5 tail better by using wider intervals.

`run_sru_timeaware.R` is the backup industrial soft-sensor intake for the
normalized MIMO SRU sequence prepared by `shared/datasets/prepare_sru.R`.
Random-split smoke (`R=10`, `B=30`) gave a small PMM win for `sru_y1_static`
(**0.953** vs LSE **0.962**) and a near tie for `sru_y1_dynamic` (PMM
**0.238**, ridge-LSE **0.237**). Chronological holdout is stricter: dynamic
`y1` has PMM **0.226** vs LSE/ridge-LSE **0.228**, while dynamic `y2` is an
almost exact PMM/LSE/ridge tie around **0.186**. SRU is therefore a useful
secondary diagnostic, not a clean flagship. Ukrainian analysis:
`../sru-c5-report.md`.

## Status & next steps

ROADMAP Phase C: C1–C3 (scaffold, core, sanity), C4 H2 benchmark, and the
GMDH-level recovery / H4 study done, including `ridge-LSE`, Huber, and L1
baselines. Ukrainian Markdown analyses: `../robust-baselines-report.md` and
`../coverage-h1-report.md`, plus C5 candidate-stage reports
`../realworld-c5-report.md`, `../domain-c5-report.md`, and
`../external-c5-fremtpl2-report.md`, `../fremtpl2-sensitivity-report.md`, and
`../fremtpl2-reserve-criterion-report.md`, and
`../fremtpl2-tail-sweep-report.md`, plus
`../gas-turbine-c5-report.md`, `../gas-turbine-timeaware-report.md`,
`../gas-turbine-arx-spike-report.md`, and
`../gas-turbine-arx-ablation-report.md`, and
`../gas-turbine-interval-coverage-report.md`, plus `../sru-c5-report.md`.
Pending: another backup external wave only if a clean flagship is still needed,
dynamic Wiener/Hammerstein systems, and a Python port. Known
limitation: on very small partitions (n ≲ 100) the PMM3 path can fire on
small-sample kurtosis bias — a candidate for a future bias-corrected gate.
