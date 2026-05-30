#!/usr/bin/env Rscript
# E1 - RUNTIME / SCALABILITY PROFILING (reviewer Q5).
#
# Question: how much wall-clock overhead does the PMM dispatch add to a full
# GMDH cascade relative to a classical LSE-GMDH cascade, and how does that
# overhead scale with the sample size n? Secondary: where does the per-candidate
# cost actually go (LSE warm-start vs bootstrap cumulant diagnostics vs PMM2
# Newton refit)?
#
# Design:
#   * KG-2-modelable regression: D=6 standardized predictors, response is the
#     degree-2 polynomial target f_true (copied from run_cascade.R, active on
#     x1..x4) plus centered skewed noise rexp(n)-1 -> gamma3 ~= 2 so the
#     dispatch genuinely routes candidates to PMM2.
#   * For n in {500,1000,2000,5000}: grow the FULL cascade (L_max=3, F=6,
#     epsilon=-Inf) twice -- force_method="LSE" (B=0) vs force_method="auto"
#     (B=200, the dispatch). Time each with system.time over 3 reps; report the
#     median wall-clock (elapsed).
#   * For n=2000 only: time the three cost components for ONE KG-2 candidate
#     (a) LSE warm-start fit, (b) one bootstrap_cumulant_diag(B=200),
#     (c) PMM2 Newton fit -- 5 reps each, median.
#
# Run from the package dir:  Rscript experiments/revision/e1_runtime.R

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

set.seed(1)

D       <- 6L
N_GRID  <- c(500L, 1000L, 2000L, 5000L)
L_MAX   <- 3L
F_TOP   <- 6L
B_PMM   <- 200L
N_REPS  <- 3L     # cascade timing reps (median wall-clock)
N_COMP  <- 5L     # per-candidate component timing: outer reps (take median)
N_INNER <- 200L   # per-candidate component timing: inner calls per rep
                  # (sub-ms resolution: time N_INNER calls, divide by N_INNER)

# Known degree-2 polynomial target (copied from run_cascade.R f_true): active on
# x1..x4; x5,x6 are decoys the cascade should ignore.
f_true <- function(X) {
  1 + 1.5 * X[, 1] - X[, 2] + 0.8 * X[, 1] * X[, 2] +
    0.6 * X[, 1]^2 - 0.4 * X[, 2]^2 + 0.5 * X[, 3] * X[, 4]
}

std <- function(x) (x - mean(x)) / stats::sd(x)

# Standardized predictors on the training matrix (per E1 spec).
make_data <- function(n, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * D), n, D)
  X <- apply(X, 2, std)                  # standardize predictors
  f <- f_true(X)
  y <- f + std(rexp(n) - 1)              # centered skewed noise (gamma3 ~ 2)
  list(X = X, y = y)
}

ctrl_lse  <- gmdh_pmm_control(B = 0L,     L_max = L_MAX, F = F_TOP,
                              epsilon = -Inf, force_method = "LSE")
ctrl_auto <- gmdh_pmm_control(B = B_PMM,  L_max = L_MAX, F = F_TOP,
                              epsilon = -Inf, force_method = "auto")

# median elapsed (wall-clock) of fitting the full cascade, over N_REPS reps
median_elapsed <- function(X, y, control) {
  ts <- numeric(N_REPS)
  for (r in seq_len(N_REPS)) {
    ts[r] <- as.numeric(system.time(gmdh_pmm(X, y, control))["elapsed"])
  }
  stats::median(ts)
}

# PMM2 share across cascade for the auto arm (sanity: confirm dispatch fires)
pmm2_share <- function(fit) {
  num <- 0; den <- 0
  for (lay in fit$layers) {
    m <- lay$methods
    den <- den + sum(m)
    if ("PMM2" %in% names(m)) num <- num + as.integer(m["PMM2"])
  }
  if (den > 0) num / den else NA_real_
}

cat(sprintf(
  "=== E1 runtime / scalability (D=%d, L_max=%d, F=%d, B_pmm=%d, reps=%d) ===\n\n",
  D, L_MAX, F_TOP, B_PMM, N_REPS))

t_global <- Sys.time()
rows <- vector("list", length(N_GRID))
for (i in seq_along(N_GRID)) {
  n <- N_GRID[i]
  dat <- make_data(n, seed = 100L + i)

  t_lse  <- median_elapsed(dat$X, dat$y, ctrl_lse)
  t_pmm  <- median_elapsed(dat$X, dat$y, ctrl_auto)

  fit_auto <- gmdh_pmm(dat$X, dat$y, ctrl_auto)   # for dispatch sanity share
  share <- pmm2_share(fit_auto)

  rows[[i]] <- data.frame(
    n           = n,
    t_lse_s     = round(t_lse, 4),
    t_pmm_s     = round(t_pmm, 4),
    overhead_x  = round(t_pmm / t_lse, 2),
    pmm2_share  = round(share, 3),
    stringsAsFactors = FALSE
  )
  cat(sprintf("  n=%5d : t_lse=%.3fs  t_pmm=%.3fs  overhead=%.2fx  PMM2share=%.2f\n",
              n, t_lse, t_pmm, t_pmm / t_lse, share))
}
scal <- do.call(rbind, rows)

cat("\nScalability table (median wall-clock over 3 reps):\n")
print(scal, row.names = FALSE)

# --- per-candidate component breakdown at n=2000 ---------------------------
cat(sprintf("\n=== per-candidate component breakdown @ n=%d (one KG-2 candidate) ===\n", 2000L))
datc <- make_data(2000L, seed = 999L)
# Use the first standardized predictor pair as the KG-2 candidate inputs.
v1 <- datc$X[, 1]; v2 <- datc$X[, 2]; y <- datc$y
d  <- kg2_design(v1, v2, y)

# Per-call median time of a thunk, at sub-ms resolution: each outer rep times
# N_INNER calls and divides by N_INNER; we take the median over N_COMP reps.
med_per_call <- function(thunk) {
  ts <- numeric(N_COMP)
  for (r in seq_len(N_COMP)) {
    el <- as.numeric(system.time(
      for (k in seq_len(N_INNER)) thunk())["elapsed"])
    ts[r] <- el / N_INNER
  }
  stats::median(ts)
}

fit_lse <- stats::lm(.kg2_formula, data = d)
eps <- stats::residuals(fit_lse)

# (a) LSE warm-start fit
t_a <- med_per_call(function() stats::lm(.kg2_formula, data = d))

# (b) one bootstrap_cumulant_diag with B=200
t_b <- med_per_call(function()
  bootstrap_cumulant_diag(eps, B = B_PMM, robust = TRUE))

# (c) PMM2 Newton fit
t_c <- med_per_call(function()
  EstemPMM::lm_pmm2(.kg2_formula, data = d, verbose = FALSE))

comp <- data.frame(
  component = c("a_lse_warmstart", "b_bootstrap_cumulant_B200", "c_pmm2_newton"),
  t_ms      = round(1000 * c(t_a, t_b, t_c), 4),
  pct_of_total = round(100 * c(t_a, t_b, t_c) / sum(t_a, t_b, t_c), 1),
  stringsAsFactors = FALSE
)
cat(sprintf("  (a) LSE warm-start            : %.4f ms\n", 1000 * t_a))
cat(sprintf("  (b) bootstrap cumulant (B=200): %.4f ms\n", 1000 * t_b))
cat(sprintf("  (c) PMM2 Newton fit           : %.4f ms\n", 1000 * t_c))
cat("\nComponent table:\n")
print(comp, row.names = FALSE)

# --- persist ---------------------------------------------------------------
outdir <- "experiments/revision"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
csv_path <- file.path(outdir, "e1_runtime.csv")

# combine into a single tidy CSV: scalability rows + component rows
scal_out <- data.frame(
  section = "scalability", n = scal$n,
  t_lse_s = scal$t_lse_s, t_pmm_s = scal$t_pmm_s,
  overhead_x = scal$overhead_x, pmm2_share = scal$pmm2_share,
  component = NA_character_, t_s = NA_real_, pct_of_total = NA_real_,
  stringsAsFactors = FALSE)
comp_out <- data.frame(
  section = "component_n2000", n = 2000L,
  t_lse_s = NA_real_, t_pmm_s = NA_real_, overhead_x = NA_real_,
  pmm2_share = NA_real_,
  component = comp$component, t_s = comp$t_ms / 1000, pct_of_total = comp$pct_of_total,
  stringsAsFactors = FALSE)
utils::write.csv(rbind(scal_out, comp_out), csv_path, row.names = FALSE)

cat(sprintf("\nSaved %s | elapsed %.1fs\n", csv_path,
            as.numeric(difftime(Sys.time(), t_global, units = "secs"))))
