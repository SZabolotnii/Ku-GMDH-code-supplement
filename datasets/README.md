# Datasets for the real-data regime screen (Paper 1A)

`external_candidates.csv` is the manifest (id, label, path, formula, source, license) used by
`experiments/run_realworld.R`. Only `concrete.csv` is shipped here (UCI Concrete Compressive
Strength, Yeh 1998; CC BY 4.0). The other benchmarks load from R packages at run time:

- airquality (ozone) — `datasets` (base R)
- Cars93           — `MASS`
- solder           — `survival`
- Credit           — `ISLR2` (a *simulated* benchmark)

Result summaries for the screen are in `experiments/results/realworld_c5_*_summary.csv`.

Note: `run_realworld.R` locates the manifest via the full monorepo layout
(`shared/datasets/external/`); in this stand-alone supplement adjust the manifest path to
`datasets/external/` (or set the path explicitly) before running the external candidates.
