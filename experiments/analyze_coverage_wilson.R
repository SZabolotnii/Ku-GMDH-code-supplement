#!/usr/bin/env Rscript
# Statistical coverage validation (Paper 1 revision, concern M6).
#
# Reads a *_raw.csv from a repeated-seed ReserveTW4 audit and, for each
# (criterion, method), tests empirical 95% prediction-interval coverage against
# the nominal 0.95 with a Wilson score interval on the pooled binomial counts
# (sum of in-interval / sum of test points across seeds). A method "passes" if
# 0.95 lies inside its Wilson 95% CI.
#
# Usage: Rscript experiments/analyze_coverage_wilson.R <raw.csv> [nominal]

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args) >= 1) args[1] else "experiments/results/volve_reservetw4_repeated_raw.csv"
nominal <- if (length(args) >= 2) as.numeric(args[2]) else 0.95
d <- utils::read.csv(path, stringsAsFactors = FALSE)

wilson <- function(x, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2); phat <- x / n
  den <- 1 + z^2 / n; ctr <- (phat + z^2 / (2 * n)) / den
  hw <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / den
  c(lower = ctr - hw, upper = ctr + hw)
}

cat(sprintf("=== Wilson coverage validity vs nominal %.2f | %s ===\n\n", nominal, basename(path)))
res <- by(d, list(d$criterion, d$method), function(g) {
  x <- sum(g$n_in); n <- sum(g$n_test); ph <- x / n
  ci <- wilson(x, n, 0.95)
  pass <- nominal >= ci["lower"] && nominal <= ci["upper"]
  verdict <- if (pass) "OK (nominal inside CI)" else if (ph < nominal) "UNDERCOVERS" else "OVERCOVERS"
  data.frame(criterion = g$criterion[1], method = g$method[1], n_pts = n,
    emp_cov = round(ph, 4), wilson_lo = round(ci["lower"], 4),
    wilson_hi = round(ci["upper"], 4), verdict = verdict, stringsAsFactors = FALSE)
})
res <- do.call(rbind, res)
res <- res[order(res$criterion, -res$emp_cov), ]
print(res, row.names = FALSE)
utils::write.csv(res, sub("_raw\\.csv$", "_wilson.csv", path), row.names = FALSE)
cat(sprintf("\nSaved %s\n", basename(sub("_raw\\.csv$", "_wilson.csv", path))))
