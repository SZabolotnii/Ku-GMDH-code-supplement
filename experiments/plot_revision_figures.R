#!/usr/bin/env Rscript
# Regenerate the revision figures (fig3 drilling NRMSE, fig4 synthetic mechanism)
# from the frozen result CSVs. Base graphics only.

resdir <- "experiments/results"
figdir_arxiv <- "../arxiv/figures"; figdir_latex <- "../latex/figures"
for (dd in c(figdir_arxiv, figdir_latex)) if (!dir.exists(dd)) dir.create(dd, recursive = TRUE)
save_both <- function(stem, draw, w = 7, h = 4) {
  for (dd in c(figdir_arxiv, figdir_latex)) {
    pdf(file.path(dd, stem), width = w, height = h); draw(); dev.off()
  }
}

methods <- c("LSE","ridge-LSE","Huber","L1","PMM")
cols <- c(LSE="#d62728", `ridge-LSE`="#ff7f0e", Huber="#2ca02c", L1="#1f77b4", PMM="#9467bd")

## ---- fig3: reserve NRMSE boxplots, Volve + FORGE (ReserveTW4) ----
v <- read.csv(file.path(resdir, "volve_reservetw4_repeated_raw.csv"))
f <- read.csv(file.path(resdir, "forge_reservetw4_repeated_raw.csv"))
v <- v[v$criterion == "ReserveTW4", ]; f <- f[f$criterion == "ReserveTW4", ]

save_both("fig3_drilling_nrmse.pdf", function() {
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  for (panel in list(list(d = v, t = "Volve 15/9-F-15 (skewed)"),
                     list(d = f, t = "Utah FORGE 58-32 (near-Gaussian)"))) {
    d <- panel$d
    bl <- lapply(methods, function(m) log10(d$nrmse_test[d$method == m]))
    names(bl) <- methods
    boxplot(bl, log = "", col = cols[methods], las = 2, ylab = "log10 reserve NRMSE",
            main = panel$t, outpch = 20, outcex = 0.5)
    abline(h = 0, lty = 3, col = "grey50")  # NRMSE = 1
  }
}, w = 9, h = 4.2)

## ---- fig4: synthetic mechanism (skewness sweep + IVH) ----
sw <- read.csv(file.path(resdir, "synthetic_var_sweep_summary.csv"))
sw <- sw[sw$criterion == "ReserveTW4", ]
agg <- aggregate(nrmse_med ~ gamma3 + method, data = sw, FUN = median)
g3 <- sort(unique(agg$gamma3))
lse <- sapply(g3, function(g) agg$nrmse_med[agg$gamma3 == g & agg$method == "LSE"])
pmm <- sapply(g3, function(g) agg$nrmse_med[agg$gamma3 == g & agg$method == "PMM"])
au <- read.csv(file.path(resdir, "synthetic_var_audit_raw.csv"))
au <- au[au$criterion == "ReserveTW4", ]

save_both("fig4_synthetic_mechanism.pdf", function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  # left: PMM advantage vs skewness
  adv <- 100 * (lse - pmm) / lse
  plot(g3, adv, type = "b", pch = 19, col = "#9467bd", lwd = 2,
       xlab = expression(gamma[3]~"(innovation skewness)"),
       ylab = "PMM advantage over LSE (% NRMSE)", main = "PMM gain scales with skewness")
  abline(h = 0, lty = 3, col = "grey50")
  # right: reserve NRMSE by method (audit), annotate IVH
  bl <- lapply(methods, function(m) au$nrmse_test[au$method == m])
  names(bl) <- methods
  boxplot(bl, col = cols[methods], las = 2, ylab = "reserve NRMSE",
          main = "Structural-break audit", ylim = c(0, 3), outpch = 20, outcex = 0.5)
  legend("topright", bty = "n", cex = 0.8,
         legend = c("66% of reserve rows", "outside training hull", "(Mahalanobis ratio 10.4x)"))
}, w = 9, h = 4.2)

cat("Wrote fig3_drilling_nrmse.pdf, fig4_synthetic_mechanism.pdf to arxiv/ and latex/ figures.\n")
