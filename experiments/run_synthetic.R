#!/usr/bin/env Rscript
# Synthetic benchmark for GMDH-PMM, Level 2 / H2 (experiment-protocol.md 3).
#
# Estimator-level Monte Carlo on a single KG-2 partial model: this isolates the
# efficiency claim (H2) from GMDH structure-selection noise, and is the main
# argument for theoretical reviewers. For each noise scenario we estimate the
# KG-2 coefficients R times with (a) forced LSE, (b) robust baselines, and
# (c) the automatic PMM dispatch, and compare empirical efficiency
#       ARE_emp    = MSE(theta_LSE) / MSE(theta_PMM)
# against theory  ARE_theory = 1 / g   (g from population cumulants).
#
# H2 pass (3.6): on PMM-favorable cells ARE_emp >= 1.1 and Wilcoxon p < 0.001;
# across all cells Pearson rho(ARE_emp, ARE_theory) > 0.85.
#
# Run from the package directory:  Rscript experiments/run_synthetic.R [R] [n] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 200L
N     <- if (length(args) >= 2) as.integer(args[2]) else 500L
B     <- if (length(args) >= 3) as.integer(args[3]) else 200L
SIGMA <- 1.0
THETA_TRUE <- c(0.5, 1.5, -1.0, 0.8, 0.6, -0.4)  # canonical KG-2 order

# --- standardized noise generators (mean 0, variance 1) --------------------
std <- function(x) (x - mean(x)) / stats::sd(x)
gens <- list(
  G        = function(n) rnorm(n),                                   # Gaussian
  exp      = function(n) rexp(n) - 1,                                # strong + skew
  lognorm  = function(n) std(rlnorm(n, 0, 0.5)),                     # A3-like skew
  uniform  = function(n) runif(n, -sqrt(3), sqrt(3)),               # B3 platykurtic
  bimodal  = function(n) std(sample(c(-1, 1), n, TRUE) * 1.3 + rnorm(n, 0, 0.45)),
  laplace  = function(n) (rexp(n) - rexp(n)) / sqrt(2),             # B1 leptokurtic (decoy)
  contam   = function(n) std(ifelse(runif(n) < 0.9, rnorm(n), rnorm(n, 0, 5)))  # C1 (decoy)
)

# population cumulants + theoretical favored method / ARE -------------------
pop_target <- function(gen, seed = 1) {
  set.seed(seed)
  sc <- sample_cumulants(gen(1e6))
  g2 <- gmdhpmm:::.g2_factor(sc$gamma3, sc$gamma4)
  g3 <- gmdhpmm:::.g3_factor(sc$gamma4, sc$gamma6)
  if (abs(sc$gamma3) > 0.3)       list(m = "PMM2", are = 1 / g2, g3 = sc$gamma3, g4 = sc$gamma4)
  else if (sc$gamma4 < -0.7)      list(m = "PMM3", are = 1 / g3, g3 = sc$gamma3, g4 = sc$gamma4)
  else                            list(m = "LSE",  are = 1.0,     g3 = sc$gamma3, g4 = sc$gamma4)
}

ctrl_pmm <- gmdh_pmm_control(B = B)
ctrl_lse <- gmdh_pmm_control(B = 0, force_method = "LSE")
ctrl_huber <- gmdh_pmm_control(B = 0, force_method = "Huber")
ctrl_l1 <- gmdh_pmm_control(B = 0, force_method = "L1")

run_cell <- function(name, gen, seed) {
  set.seed(seed)
  err_lse <- numeric(R_REP); err_huber <- numeric(R_REP)
  err_l1 <- numeric(R_REP); err_pmm <- numeric(R_REP)
  methods <- character(R_REP)
  for (r in seq_len(R_REP)) {
    v1 <- rnorm(N); v2 <- rnorm(N)
    y <- kg2_predict(THETA_TRUE, v1, v2) + SIGMA * gen(N)
    th_lse <- inner_estimate(v1, v2, y, ctrl_lse)$theta
    th_huber <- inner_estimate(v1, v2, y, ctrl_huber)$theta
    th_l1 <- inner_estimate(v1, v2, y, ctrl_l1)$theta
    est <- inner_estimate(v1, v2, y, ctrl_pmm)
    err_lse[r] <- sum((th_lse - THETA_TRUE)^2)
    err_huber[r] <- sum((th_huber - THETA_TRUE)^2)
    err_l1[r] <- sum((th_l1 - THETA_TRUE)^2)
    err_pmm[r] <- sum((est$theta - THETA_TRUE)^2)
    methods[r] <- est$method
  }
  tgt <- pop_target(gen)
  are_emp <- mean(err_lse) / mean(err_pmm)
  wp <- suppressWarnings(stats::wilcox.test(err_lse, err_pmm, paired = TRUE)$p.value)
  data.frame(
    scenario   = name,
    gamma3     = round(tgt$g3, 3),
    gamma4     = round(tgt$g4, 3),
    favored    = tgt$m,
    dispatch   = names(sort(table(methods), decreasing = TRUE))[1],
    disp_frac  = round(mean(methods == tgt$m), 2),
    ARE_theory = round(tgt$are, 3),
    ARE_emp    = round(are_emp, 3),
    ARE_Huber  = round(mean(err_lse) / mean(err_huber), 3),
    ARE_L1     = round(mean(err_lse) / mean(err_l1), 3),
    PMM_vs_Huber = round(mean(err_huber) / mean(err_pmm), 3),
    PMM_vs_L1    = round(mean(err_l1) / mean(err_pmm), 3),
    var_red_pct = round((1 - mean(err_pmm) / mean(err_lse)) * 100, 1),
    wilcox_p   = signif(wp, 3),
    stringsAsFactors = FALSE
  )
}

cat(sprintf("=== GMDH-PMM synthetic H2 benchmark (R=%d, n=%d, B=%d) ===\n\n",
            R_REP, N, B))
t0 <- Sys.time()
res <- do.call(rbind, Map(run_cell, names(gens), gens, seq_along(gens) + 1000L))
cat("Per-scenario results:\n")
print(res, row.names = FALSE)

# --- H2 verdict ------------------------------------------------------------
rho <- stats::cor(res$ARE_emp, res$ARE_theory)
fav <- res[res$favored != "LSE", ]
h2_are <- all(fav$ARE_emp >= 1.1)
h2_sig <- all(fav$wilcox_p < 0.001)
h2_rho <- rho > 0.85

cat(sprintf("\nPearson rho(ARE_emp, ARE_theory) = %.3f  (pass > 0.85: %s)\n", rho, h2_rho))
cat(sprintf("PMM-favorable cells: ARE_emp >= 1.1 : %s | Wilcoxon p < 1e-3 : %s\n",
            h2_are, h2_sig))
cat(sprintf(">>> H2 %s\n", if (h2_are && h2_sig && h2_rho) "PASS" else "PARTIAL/FAIL"))

# --- persist ---------------------------------------------------------------
outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(res, file.path(outdir, "synthetic_h2.csv"), row.names = FALSE)

# ARE_emp vs ARE_theory scatter (H2 visual): points should track the y = x line.
png(file.path(outdir, "are_emp_vs_theory.png"), width = 720, height = 640, res = 110)
lim <- range(0.8, res$ARE_theory, res$ARE_emp)
plot(res$ARE_theory, res$ARE_emp, xlim = lim, ylim = lim, pch = 19,
     col = ifelse(res$favored == "LSE", "grey50", "firebrick"),
     xlab = "ARE_theory = 1 / g", ylab = "ARE_emp = MSE(LSE)/MSE(PMM)",
     main = sprintf("GMDH-PMM efficiency vs theory (rho = %.3f)", rho))
abline(0, 1, lty = 2, col = "grey40")
text(res$ARE_theory, res$ARE_emp, res$scenario, pos = 4, cex = 0.7)
legend("topleft", c("PMM-favorable", "LSE (decoy)"), pch = 19,
       col = c("firebrick", "grey50"), bty = "n", cex = 0.8)
invisible(dev.off())

cat(sprintf("\nSaved %s (+ are_emp_vs_theory.png) | elapsed %.1fs\n",
            file.path(outdir, "synthetic_h2.csv"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
