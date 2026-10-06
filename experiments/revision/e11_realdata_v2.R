#!/usr/bin/env Rscript
# E11 -- Real-data regime screen, corrected protocol (R1 prep, 2026-10-06).
#
# Fixes three defects of the submitted Tab. realdata:
#   1. the cells were MEANS over splits under a "median" caption;
#   2. each estimator arm grew its cascade on a different inner train/val split
#      (seed + arm index); the _v2 runs share one inner seed per split;
#   3. the 30 repeated 75/25 splits of one dataset are not independent, and the
#      5 datasets x 4 competitors = 20 comparisons were not corrected.
# For every dataset and competitor X it reports the median NRMSE, the paired
# median difference PMM - X, the Wilcoxon p, the Nadeau-Bengio corrected
# resampled t-test p (variance inflated by 1/J + n_test/n_train), and Holm
# adjustment of both over the 20 comparisons.
#
# Inputs : experiments/results/realworld_c5_<id>_v2_raw.csv  (run_realworld.R <id> 30 200 _v2)
# Output : experiments/revision/e11_realdata_v2.csv
# Run from the package directory: Rscript experiments/revision/e11_realdata_v2.R

IDS <- c(Credit = "islr_credit_balance", Concrete = "concrete",
         airquality = "airquality_ozone", Cars93 = "mass_cars93_mpg",
         solder = "survival_solder_skips")
COMPETITORS <- c("LSE", "ridge-LSE", "Huber", "L1")
TEST_FRAC <- 0.25

rows <- list()
for (nm in names(IDS)) {
  d <- utils::read.csv(sprintf("experiments/results/realworld_c5_%s_v2_raw.csv", IDS[[nm]]),
                       stringsAsFactors = FALSE)
  d <- d[d$ok, ]
  p <- d[d$method == "PMM", ]
  med <- tapply(d$nrmse, d$method, stats::median)
  for (m in COMPETITORS) {
    x <- merge(p[, c("rep", "nrmse")], d[d$method == m, c("rep", "nrmse")], by = "rep")
    dif <- x$nrmse.x - x$nrmse.y; J <- length(dif)
    se_nb <- sqrt((1 / J + TEST_FRAC / (1 - TEST_FRAC)) * stats::var(dif))
    t_nb <- if (se_nb > 0) mean(dif) / se_nb else 0
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, competitor = m, splits = J,
      pmm2_share = round(mean(p$pmm2_share, na.rm = TRUE), 2),
      median_pmm = round(med[["PMM"]], 3), median_competitor = round(med[[m]], 3),
      median_diff = round(stats::median(dif), 4), pmm_wins = sum(dif < 0),
      wilcox_p = suppressWarnings(stats::wilcox.test(x$nrmse.x, x$nrmse.y, paired = TRUE)$p.value),
      nb_p = 2 * stats::pt(-abs(t_nb), df = J - 1),
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
res$wilcox_p_holm <- stats::p.adjust(res$wilcox_p, "holm")
res$nb_p_holm <- stats::p.adjust(res$nb_p, "holm")
for (k in c("wilcox_p", "nb_p", "wilcox_p_holm", "nb_p_holm")) res[[k]] <- signif(res[[k]], 3)
options(width = 220)
print(res, row.names = FALSE)
utils::write.csv(res, "experiments/revision/e11_realdata_v2.csv", row.names = FALSE)
cat("\nSaved experiments/revision/e11_realdata_v2.csv\n")
