#!/usr/bin/env Rscript
# E10 -- Seed-clustered coverage for the synthetic VAR audit (R1 prep, 2026-10-06).
#
# The manuscript's coverage paragraph pools in-interval indicators over all test
# points of all 30 seeds into one binomial and puts a Wilson interval on it.
# Test points within a seed share one fitted cascade and one calibration
# quantile, so they are not independent Bernoulli draws: the pooled Wilson
# interval is anti-conservative by the design effect. Here the seed is the unit:
#   - per-seed coverage c_s; mean over seeds; 95% CI by seed bootstrap and t;
#   - design effect = Var_seed(c_s) / (p(1-p)/n_test)  (>1 => pooling too narrow);
#   - paired per-seed coverage differences PMM - X (Wilcoxon), so "closest to
#     nominal" is a tested statement, not an ordering of point estimates.
#
# Input : experiments/results/synthetic_var_audit_raw.csv (ReserveTW4 rows)
# Output: experiments/revision/e10_coverage_cluster.csv
# Run from the package directory: Rscript experiments/revision/e10_coverage_cluster.R

set.seed(20261006L)
d <- utils::read.csv("experiments/results/synthetic_var_audit_raw.csv", stringsAsFactors = FALSE)
d <- d[d$criterion == "ReserveTW4", ]
N_TEST <- 2500L   # reserve window = post-break half of N = 5000 (sec:varbreak)
NOMINAL <- 0.95

wilson <- function(p, n, z = 1.96) {
  den <- 1 + z^2 / n; ctr <- (p + z^2 / (2 * n)) / den
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(ctr - hw, ctr + hw)
}

rows <- lapply(split(d, d$method), function(g) {
  c_s <- g$coverage; S <- length(c_s); p <- mean(c_s)
  bt <- replicate(5000, mean(sample(c_s, S, TRUE)))
  tci <- p + c(-1, 1) * stats::qt(0.975, S - 1) * stats::sd(c_s) / sqrt(S)
  w <- wilson(p, S * N_TEST)
  data.frame(method = g$method[1], seeds = S, coverage_mean = round(p, 4),
             coverage_median = round(stats::median(c_s), 4),
             seed_boot_lo = round(stats::quantile(bt, 0.025, names = FALSE), 4),
             seed_boot_hi = round(stats::quantile(bt, 0.975, names = FALSE), 4),
             seed_t_lo = round(tci[1], 4), seed_t_hi = round(tci[2], 4),
             pooled_wilson_lo = round(w[1], 4), pooled_wilson_hi = round(w[2], 4),
             design_effect = round(stats::var(c_s) / (p * (1 - p) / N_TEST), 1),
             nominal_in_seed_ci = NOMINAL >= tci[1] && NOMINAL <= tci[2],
             stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)
res <- res[order(-res$coverage_mean), ]

pmm <- d[d$method == "PMM", c("seed", "coverage")]
res$paired_diff_pmm_minus <- NA_real_; res$paired_p <- NA_real_
for (m in setdiff(res$method, "PMM")) {
  x <- merge(pmm, d[d$method == m, c("seed", "coverage")], by = "seed")
  i <- res$method == m
  res$paired_diff_pmm_minus[i] <- round(stats::median(x$coverage.x - x$coverage.y), 4)
  res$paired_p[i] <- signif(suppressWarnings(
    stats::wilcox.test(x$coverage.x, x$coverage.y, paired = TRUE)$p.value), 3)
}
options(width = 200)
print(res, row.names = FALSE)
utils::write.csv(res, "experiments/revision/e10_coverage_cluster.csv", row.names = FALSE)
cat("\nSaved experiments/revision/e10_coverage_cluster.csv\n")
