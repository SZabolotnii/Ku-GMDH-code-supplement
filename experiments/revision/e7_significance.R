#!/usr/bin/env Rscript
# E7 -- Significance tests / CIs for method differences.
#
# Reviewer request: report tests beyond medians/IQRs for the headline
# PMM-vs-baseline comparisons. For each dataset and each baseline X we compute,
# pairing PMM and X by the same unit (seed or rep):
#   - the paired MEDIAN of the per-pair NRMSE difference (PMM - X);
#     negative => PMM is better (lower error).
#   - a 95% bootstrap CI for that median difference (2000 resamples over pairs).
#   - the paired Wilcoxon signed-rank p-value (exact when possible).
#
# Inputs (existing raw CSVs, no re-running of experiments):
#   experiments/results/synthetic_var_audit_raw.csv      (criterion==ReserveTW4)
#   experiments/results/realworld_c5_islr_credit_balance_raw.csv
# A cascade raw CSV does NOT exist (cascade_recovery.csv is pre-aggregated
# means/sds with no per-seed values), so no paired test is possible there.
#
# Output:
#   experiments/revision/e7_significance.csv

suppressMessages({
  library(stats)
})

set.seed(20260530L)

RESULTS_DIR <- "experiments/results"
OUT_CSV     <- "experiments/revision/e7_significance.csv"
B_BOOT      <- 2000L
BASELINES   <- c("LSE", "ridge-LSE", "Huber", "L1")

## --- helpers ---------------------------------------------------------------

# Bootstrap CI for the median of a paired difference vector d (= value_PMM - value_X).
boot_median_ci <- function(d, B = B_BOOT, conf = 0.95) {
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 2L) return(c(low = NA_real_, high = NA_real_))
  stats <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    stats[b] <- median(d[idx])
  }
  a <- (1 - conf) / 2
  qs <- quantile(stats, probs = c(a, 1 - a), names = FALSE, type = 7)
  c(low = qs[1], high = qs[2])
}

# One PMM-vs-X paired comparison given a wide table with a 'unit' id column,
# a PMM numeric column and an X numeric column.
paired_compare <- function(unit, pmm, x) {
  ok <- is.finite(pmm) & is.finite(x)
  pmm <- pmm[ok]; x <- x[ok]
  d <- pmm - x                       # negative => PMM lower NRMSE => PMM better
  n <- length(d)
  med <- median(d)
  ci  <- boot_median_ci(d)
  # paired Wilcoxon signed-rank; exact if no ties/zeros and n small enough
  wp <- tryCatch(
    suppressWarnings(wilcox.test(pmm, x, paired = TRUE)$p.value),
    error = function(e) NA_real_
  )
  list(median_diff = med, ci_low = ci["low"], ci_high = ci["high"],
       wilcox_p = wp, n_pairs = n)
}

# Run all PMM-vs-baseline comparisons for one dataset.
# df: long data with columns unit_col, method_col, value_col.
run_dataset <- function(df, dataset, unit_col, method_col, value_col,
                        baselines = BASELINES) {
  # wide: one row per unit, one column per method
  w <- reshape(
    df[, c(unit_col, method_col, value_col)],
    idvar = unit_col, timevar = method_col, direction = "wide"
  )
  colnames(w) <- sub(paste0("^", value_col, "\\."), "", colnames(w))
  stopifnot("PMM" %in% colnames(w))
  rows <- list()
  for (bx in baselines) {
    if (!bx %in% colnames(w)) {
      message(sprintf("  [%s] baseline '%s' not present -- skipped", dataset, bx))
      next
    }
    r <- paired_compare(w[[unit_col]], w[["PMM"]], w[[bx]])
    rows[[length(rows) + 1L]] <- data.frame(
      dataset     = dataset,
      comparison  = sprintf("PMM vs %s", bx),
      median_diff = r$median_diff,
      ci_low      = unname(r$ci_low),
      ci_high     = unname(r$ci_high),
      wilcox_p    = r$wilcox_p,
      n_pairs     = r$n_pairs,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

## --- (1) Synthetic VAR audit, criterion = ReserveTW4 -----------------------

va_path <- file.path(RESULTS_DIR, "synthetic_var_audit_raw.csv")
va <- read.csv(va_path, stringsAsFactors = FALSE)
va <- va[va$criterion == "ReserveTW4", ]
cat(sprintf("[VAR audit] criterion=ReserveTW4: %d rows, %d seeds, methods: %s\n",
            nrow(va), length(unique(va$seed)), paste(unique(va$method), collapse = ", ")))
res_var <- run_dataset(va, "synthetic_var_audit (ReserveTW4)",
                       unit_col = "seed", method_col = "method",
                       value_col = "nrmse_test")

## --- (2) ISLR Credit balance, paired by rep --------------------------------

cr_path <- file.path(RESULTS_DIR, "realworld_c5_islr_credit_balance_raw.csv")
cr <- read.csv(cr_path, stringsAsFactors = FALSE)
cat(sprintf("[Credit balance] %d rows, %d reps, methods: %s\n",
            nrow(cr), length(unique(cr$rep)), paste(unique(cr$method), collapse = ", ")))
res_cr <- run_dataset(cr, "islr_credit_balance",
                      unit_col = "rep", method_col = "method",
                      value_col = "nrmse")

## --- (3) Cascade -----------------------------------------------------------
# cascade_recovery.csv is pre-aggregated (rmse_mean/rmse_sd by noise/arm/depth);
# it contains NO per-seed raw values, so a paired test is not possible.
cat("[Cascade] no per-seed raw CSV exists (cascade_recovery.csv is aggregated) -- skipped.\n")

## --- combine + write -------------------------------------------------------

res <- rbind(res_var, res_cr)
# round for display readability but keep full precision in memory before write
res_disp <- res
num_cols <- c("median_diff", "ci_low", "ci_high")
res_disp[num_cols] <- lapply(res_disp[num_cols], function(z) round(z, 5))
res_disp$wilcox_p  <- signif(res_disp$wilcox_p, 4)

write.csv(res, OUT_CSV, row.names = FALSE)

cat("\n================ E7 SIGNIFICANCE TABLE ================\n")
cat("(median_diff = median(PMM - X) NRMSE; negative => PMM better)\n\n")
print(res_disp, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", normalizePath(OUT_CSV, mustWork = FALSE)))
