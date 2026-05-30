#!/usr/bin/env Rscript
# Generate submission-facing figures for the GMDH-PMM manuscript.
#
# The script uses only saved CSV artifacts, so it is safe to run after the
# computational experiments have been completed.

find_root <- function() {
  cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "ROADMAP.md")) &&
        dir.exists(file.path(cur, "paper-1-gmdh-pmm"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Cannot locate project root.")
    cur <- parent
  }
}

root <- find_root()
res_dir <- file.path(root, "paper-1-gmdh-pmm", "code", "experiments", "results")
fig_dir <- file.path(root, "paper-1-gmdh-pmm", "latex", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name) {
  path <- file.path(res_dir, name)
  if (!file.exists(path)) stop("Missing result CSV: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE)
}

method_cols <- c(PMM = "#1b9e77", LSE = "#d95f02", `ridge-LSE-1e-1` = "#66a61e",
                 Huber = "#7570b3", L1 = "#e7298a")

synthetic <- read_result("synthetic_h2.csv")
pdf(file.path(fig_dir, "fig1_synthetic_are.pdf"), width = 7.2, height = 4.6,
    family = "Helvetica")
op <- par(mar = c(7.2, 4.5, 1.2, 1.0))
sc <- synthetic$scenario
x <- seq_along(sc)
plot(x, synthetic$ARE_emp, type = "b", pch = 19, lwd = 2, col = "#1b9e77",
     xaxt = "n", ylim = c(0, max(synthetic$ARE_theory, synthetic$ARE_emp) * 1.08),
     xlab = "", ylab = "ARE = MSE(LSE) / MSE(PMM)")
lines(x, synthetic$ARE_theory, type = "b", pch = 17, lwd = 2, col = "#386cb0")
abline(h = 1, lty = 3, col = "grey45")
axis(1, at = x, labels = sc, las = 2, cex.axis = 0.85)
legend("topleft", c("empirical ARE", "theoretical ARE"), pch = c(19, 17),
       lwd = 2, col = c("#1b9e77", "#386cb0"), bty = "n")
par(op)
dev.off()

coverage <- read_result("gas_turbine_interval_coverage_co_summary.csv")
coverage <- coverage[coverage$criterion_label == "ReserveTW4" &
                       coverage$interval_type == "abs_residual_95", ]
coverage <- coverage[match(c("PMM", "LSE", "Huber", "L1"), coverage$method), ]
pdf(file.path(fig_dir, "fig2_coverage_width_tradeoff.pdf"), width = 6.4, height = 4.7,
    family = "Helvetica")
op <- par(mar = c(4.4, 4.5, 1.2, 1.0))
plot(coverage$test_mean_width, coverage$test_coverage,
     pch = 19, cex = 1.5, col = method_cols[coverage$method],
     xlab = "Mean prediction-interval width",
     ylab = "Test coverage",
     xlim = range(coverage$test_mean_width) + c(-0.15, 0.18),
     ylim = range(coverage$test_coverage, 0.95) + c(-0.012, 0.008))
abline(h = 0.95, lty = 3, col = "grey45")
text(coverage$test_mean_width, coverage$test_coverage,
     labels = paste0(coverage$method, "\nNRMSE=", sprintf("%.3f", coverage$nrmse)),
     pos = c(4, 4, 2, 2), cex = 0.78)
par(op)
dev.off()

audit <- read_result("gas_turbine_calibration_audit_co_summary.csv")
audit <- audit[audit$interval_type == "abs_residual_95", ]
audit <- audit[audit$criterion_label %in% c("MSE", "ReserveTW4"), ]
audit <- audit[audit$method %in% names(method_cols), ]
audit$method <- factor(audit$method, levels = names(method_cols))
audit <- audit[order(audit$criterion_label, audit$method), ]
pdf(file.path(fig_dir, "fig3_nested_audit_nrmse.pdf"), width = 7.0, height = 4.8,
    family = "Helvetica")
op <- par(mar = c(5.2, 4.7, 1.2, 1.0))
mat <- tapply(audit$nrmse, list(audit$method, audit$criterion_label), identity)
bp <- barplot(mat, beside = TRUE, log = "y", col = method_cols[rownames(mat)],
              ylab = "NRMSE, log scale", ylim = c(0.45, 30),
              names.arg = colnames(mat), las = 1)
legend("topleft", rownames(mat), fill = method_cols[rownames(mat)], bty = "n",
       ncol = 2, cex = 0.85)
text(bp, mat, labels = sprintf("%.3f", mat), pos = 3, cex = 0.72)
par(op)
dev.off()

reserve <- audit[audit$criterion_label == "ReserveTW4", ]
reserve <- reserve[order(reserve$method), ]
pdf(file.path(fig_dir, "fig4_reserve_stability_diagnostics.pdf"), width = 7.0, height = 4.8,
    family = "Helvetica")
op <- par(mfrow = c(1, 2), mar = c(6.2, 4.4, 1.2, 0.8))
barplot(setNames(reserve$nrmse, reserve$method), log = "y",
        col = method_cols[as.character(reserve$method)], las = 2,
        ylab = "Test NRMSE, log scale")
barplot(setNames(reserve$max_log10_condition, reserve$method),
        col = method_cols[as.character(reserve$method)], las = 2,
        ylab = expression(max~log[10]~kappa))
par(op)
dev.off()

cat("Generated manuscript figures in ", fig_dir, "\n", sep = "")
