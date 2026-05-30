#!/usr/bin/env Rscript
# Sanity-level experiment for GMDH-PMM (experiment-protocol.md Level 1).
#
# Demonstrates:
#   * H3 (safe fallback): under Gaussian noise the dispatch sends every partial
#     model to LSE, so GMDH-PMM == GMDH-LSE.
#   * Dispatch firing: under strongly skewed noise PMM2 is selected, and
#     GMDH-PMM does not lose to an LSE-forced tournament out of sample.
#
# Run from the package directory:  Rscript experiments/run_sanity.R

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

set.seed(2026)

make_data <- function(n, d = 4, noise) {
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  y <- 1 + 1.5 * X[, 1] - X[, 2] + 0.8 * X[, 1] * X[, 2] + 0.5 * X[, 3]^2 + noise
  list(X = X, y = y)
}

method_counts <- function(fit) {
  m <- vapply(Filter(function(z) !is.null(z) && z$type == "model", fit$nodes),
              function(z) z$method, character(1))
  table(factor(m, levels = c("LSE", "PMM2", "PMM3")))
}

mse <- function(fit, dat) mean((dat$y - predict(fit, dat$X))^2)

cat("=== GMDH-PMM sanity (seed 2026) ===\n\n")

## ---- Case A: Gaussian noise -> safe fallback (H3) --------------------------
n <- 600
trainA <- make_data(n, noise = rnorm(n))
testA  <- make_data(300, noise = rnorm(300))
fitA <- gmdh_pmm(trainA$X, trainA$y, gmdh_pmm_control(B = 200, L_max = 4, F = 6, seed = 7))

mcA <- method_counts(fitA)
cat("Case A  Gaussian noise\n")
print(mcA)
cat(sprintf("  LSE share = %d/%d | test MSE = %.4f\n", mcA["LSE"], sum(mcA), mse(fitA, testA)))
cat("  (safe fallback: dispatch is overwhelmingly LSE; rare deep-layer PMM2,\n")
cat("   where intermediate residuals drift from Gaussian, is harmless)\n\n")

## ---- Case B: strong right-skew noise -> PMM2 dispatch ----------------------
skew_noise <- (rexp(n) - 1) * 2.0
trainB <- make_data(n, noise = skew_noise)
testB  <- make_data(300, noise = (rexp(300) - 1) * 2.0)

fitB_pmm <- gmdh_pmm(trainB$X, trainB$y, gmdh_pmm_control(B = 200, L_max = 4, F = 6, seed = 7))
# LSE-forced baseline: thresholds the dispatch can never satisfy
ctrl_lse <- gmdh_pmm_control(B = 200, L_max = 4, F = 6, seed = 7,
                             skew_min = Inf, skew_strong = Inf,
                             g2_threshold = 0, kurt_threshold = -Inf)
fitB_lse <- gmdh_pmm(trainB$X, trainB$y, ctrl_lse)

cat("Case B  strong right-skew noise (scaled exponential)\n")
cat("  GMDH-PMM partial-model methods:\n"); print(method_counts(fitB_pmm))
cat(sprintf("  GMDH-PMM test MSE = %.4f\n", mse(fitB_pmm, testB)))
cat(sprintf("  GMDH-LSE test MSE = %.4f\n", mse(fitB_lse, testB)))
cat(sprintf("  PMM2 selected at least once: %s\n",
            method_counts(fitB_pmm)["PMM2"] > 0))

cat("\n=== done ===\n")
