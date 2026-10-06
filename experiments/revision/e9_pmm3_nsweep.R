#!/usr/bin/env Rscript
# E9 -- n-sweep of the PMM3 theory-vs-empirical ARE gap (R1 prep, 2026-10-06).
#
# Tab. synthetic reports, at n=2000, ARE_theory 3.334 vs ARE_emp 2.697 (uniform)
# and 6.122 vs 3.100 (bimodal). The manuscript attributes the gap to slow
# convergence of the sixth-cumulant term. If that is right, ARE_emp must climb
# toward 1/g3 as n grows. This script measures it on the same single KG-2
# partial model as run_synthetic.R, with PMM3 FORCED (no dispatch noise):
#   ARE_emp(n) = MSE(theta_LSE) / MSE(theta_PMM3),  paired by replication.
#
# Run from the package directory:
#   Rscript experiments/revision/e9_pmm3_nsweep.R [R] [n_grid] [regimes] [out_tag]
# Output: experiments/revision/e9_pmm3_nsweep[out_tag].csv

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 400L
N_GRID <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else
  c(250L, 500L, 1000L, 2000L, 5000L, 10000L, 20000L)
THETA_TRUE <- c(0.5, 1.5, -1.0, 0.8, 0.6, -0.4)   # as in run_synthetic.R
B_BOOT <- 2000L

std <- function(x) (x - mean(x)) / stats::sd(x)
# `bimodal` is the generator of run_synthetic.R / Tab. synthetic. Its std() centres
# the SAMPLE, so the noise has sample mean exactly 0: it is orthogonal to the
# intercept column and hands LSE an intercept error far below its sampling
# variance. `bimodal_pop` standardizes with the POPULATION moments
# (mean 0, var 1.3^2 + 0.45^2), which is what the ARE theory assumes.
gens <- list(
  uniform = function(n) runif(n, -sqrt(3), sqrt(3)),
  bimodal = function(n) std(sample(c(-1, 1), n, TRUE) * 1.3 + rnorm(n, 0, 0.45)),
  bimodal_pop = function(n) (sample(c(-1, 1), n, TRUE) * 1.3 + rnorm(n, 0, 0.45)) /
    sqrt(1.3^2 + 0.45^2)
)
if (length(args) >= 3) gens <- gens[strsplit(args[3], ",")[[1]]]

ctrl_lse  <- gmdh_pmm_control(B = 0, force_method = "LSE")
ctrl_pmm3 <- gmdh_pmm_control(B = 0, force_method = "PMM3")

are_theory <- function(gen) {
  set.seed(1)
  sc <- sample_cumulants(gen(1e6))
  1 / gmdhpmm:::.g3_factor(sc$gamma4, sc$gamma6)
}

rows <- list()
t0 <- Sys.time()
for (nm in names(gens)) {
  th <- are_theory(gens[[nm]])
  for (n in N_GRID) {
    set.seed(9000L + n + match(nm, names(gens)))
    e_lse <- e_pmm <- i_lse <- i_pmm <- numeric(R_REP); fell_back <- 0L
    for (r in seq_len(R_REP)) {
      v1 <- rnorm(n); v2 <- rnorm(n)
      y <- kg2_predict(THETA_TRUE, v1, v2) + gens[[nm]](n)
      d_lse <- (inner_estimate(v1, v2, y, ctrl_lse)$theta - THETA_TRUE)^2
      est <- inner_estimate(v1, v2, y, ctrl_pmm3)
      if (!identical(est$method, "PMM3")) fell_back <- fell_back + 1L
      d_pmm <- (est$theta - THETA_TRUE)^2
      e_lse[r] <- sum(d_lse); e_pmm[r] <- sum(d_pmm)
      i_lse[r] <- d_lse[1]; i_pmm[r] <- d_pmm[1]   # canonical order: intercept first
    }
    are <- mean(e_lse) / mean(e_pmm)
    bt <- replicate(B_BOOT, { i <- sample.int(R_REP, R_REP, TRUE); mean(e_lse[i]) / mean(e_pmm[i]) })
    rows[[length(rows) + 1]] <- data.frame(
      regime = nm, n = n, reps = R_REP, ARE_theory = round(th, 3),
      ARE_emp = round(are, 3),
      ARE_lo = round(stats::quantile(bt, 0.025, names = FALSE), 3),
      ARE_hi = round(stats::quantile(bt, 0.975, names = FALSE), 3),
      ratio_emp_theory = round(are / th, 3),
      ARE_intercept = round(mean(i_lse) / mean(i_pmm), 3),
      ARE_non_intercept = round(mean(e_lse - i_lse) / mean(e_pmm - i_pmm), 3),
      pmm3_fallbacks = fell_back)
    z <- rows[[length(rows)]]
    cat(sprintf("%-11s n=%6d  ARE_emp %.3f [%.3f, %.3f]  theory %.3f (%.0f%%)  intercept %.2f  rest %.2f  fallbacks %d\n",
                nm, n, are, z$ARE_lo, z$ARE_hi, th, 100 * are / th,
                z$ARE_intercept, z$ARE_non_intercept, fell_back))
  }
}
res <- do.call(rbind, rows)
out <- sprintf("experiments/revision/e9_pmm3_nsweep%s.csv", if (length(args) >= 4) args[4] else "")
utils::write.csv(res, out, row.names = FALSE)
cat(sprintf("\nSaved %s | elapsed %.1f min\n", out,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
