# Datasets for the real-data regime screen

`external_candidates.csv` is the manifest (id, label, path, formula, source, license) used by
`experiments/run_realworld.R`. Only `concrete.csv` is shipped here (UCI Concrete Compressive
Strength, Yeh 1998; CC BY 4.0). The other benchmarks load from R packages at run time:

- airquality (ozone) — `datasets` (base R)
- Cars93           — `MASS`
- solder           — `survival`
- Credit           — `ISLR2` (a *simulated* benchmark)

Result summaries for the screen are in `experiments/results/realworld_c5_*_summary.csv`.

In this standalone supplement, `run_realworld.R` first uses
`datasets/external/external_candidates.csv` and resolves shipped public files from
`datasets/external/`. The unshipped external candidates remain in the manifest for
provenance, but their archived results should be treated as precomputed unless the
corresponding prepared raw files are added locally.
