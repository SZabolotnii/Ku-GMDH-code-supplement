#!/usr/bin/env Rscript
# EXPERIMENT E3 -- ABLATIONS (reviewer concern: which components drive the gains).
#
# Estimator-level Monte Carlo on a single KG-2 partial model (same design as
# experiments/run_synthetic.R), so the ablation isolates the inner estimator and
# the *dispatch* from GMDH structure-selection noise. For each of the 7 residual
# regimes we fit the canonical KG-2 coefficient vector R times and compare arms:
#
#   (i)   FULL   -- automatic dispatch with bootstrap-stabilized cumulants (B=200)
#   (ii)  B0     -- dispatch with NO bootstrap (plug-in sample cumulants + their
#                   Gaussian analytic SEs se3=sqrt(6/n), se4=sqrt(24/n))
#   (iii) PMM2   -- always force PMM2 (no dispatch)
#   (iv)  PMM3   -- always force PMM3 (no dispatch)
#   (v)   LSE    -- ordinary least squares baseline
#
# Headline metric per (regime, arm): empirical ARE = MSE(LSE) / MSE(arm).
#   ARE > 1  => arm beats LSE on coefficient MSE; ARE < 1 => arm worse than LSE.
# We also record, for the two dispatch arms, WHICH method was selected (the modal
# choice and the activation share of each method) -- this is the heart of the two
# claims:
#   (a) the bootstrap SUPPRESSES false PMM activation on the heavy-tailed decoys
#       (contaminated, Laplace) relative to B=0;
#   (b) always-forcing PMM (no dispatch) HURTS on Gaussian + decoy regimes versus
#       the dispatch.
#
# Usage: Rscript experiments/revision/e3_ablation.R [R] [n] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args  <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 50L
N     <- if (length(args) >= 2) as.integer(args[2]) else 2000L
B     <- if (length(args) >= 3) as.integer(args[3]) else 200L
SIGMA <- 1.0
THETA_TRUE <- c(0.5, 1.5, -1.0, 0.8, 0.6, -0.4)   # canonical KG-2 order

# --- standardized residual regimes (mean 0, variance 1) --------------------
# Same 7 regimes as run_synthetic.R: 2 PMM2-favorable (skew), 1 PMM3-favorable
# (platykurtic uniform), plus Gaussian and two HEAVY-TAILED DECOYS where PMM
# must NOT activate (laplace = leptokurtic symmetric, contam = 90/10 mixture).
std <- function(x) (x - mean(x)) / stats::sd(x)
gens <- list(
  Gaussian = function(n) rnorm(n),                                   # null
  exp      = function(n) rexp(n) - 1,                                # strong + skew (PMM2)
  lognorm  = function(n) std(rlnorm(n, 0, 0.5)),                     # A3-like skew (PMM2)
  uniform  = function(n) runif(n, -sqrt(3), sqrt(3)),               # platykurtic (PMM3)
  bimodal  = function(n) std(sample(c(-1, 1), n, TRUE) * 1.3 + rnorm(n, 0, 0.45)),
  laplace  = function(n) (rexp(n) - rexp(n)) / sqrt(2),             # leptokurtic DECOY
  contam   = function(n) std(ifelse(runif(n) < 0.9, rnorm(n), rnorm(n, 0, 5)))  # mixture DECOY
)

# population favored method (the "oracle" dispatch target) from 1e6 draws -----
pop_target <- function(gen, seed = 1) {
  set.seed(seed)
  sc <- sample_cumulants(gen(1e6))
  g2 <- gmdhpmm:::.g2_factor(sc$gamma3, sc$gamma4)
  g3 <- gmdhpmm:::.g3_factor(sc$gamma4, sc$gamma6)
  if (abs(sc$gamma3) > 0.3)  list(m = "PMM2", are = 1 / g2, g3 = sc$gamma3, g4 = sc$gamma4)
  else if (sc$gamma4 < -0.7) list(m = "PMM3", are = 1 / g3, g3 = sc$gamma3, g4 = sc$gamma4)
  else                       list(m = "LSE",  are = 1.0,    g3 = sc$gamma3, g4 = sc$gamma4)
}

# --- the five arms ----------------------------------------------------------
# Bootstrap arms get a per-replicate boot_seed for reproducibility.
ctrl_full <- function(seed) gmdh_pmm_control(B = B, force_method = "auto", boot_seed = seed)
ctrl_b0   <-                gmdh_pmm_control(B = 0, force_method = "auto")
ctrl_pmm2 <-                gmdh_pmm_control(B = 0, force_method = "PMM2")
ctrl_pmm3 <-                gmdh_pmm_control(B = 0, force_method = "PMM3")
ctrl_lse  <-                gmdh_pmm_control(B = 0, force_method = "LSE")

ARMS <- c("FULL", "B0", "PMM2", "PMM3", "LSE")

run_cell <- function(name, gen, seed) {
  set.seed(seed)
  # squared-coefficient-error accumulators + realized methods per arm
  err <- lapply(ARMS, function(a) numeric(R_REP)); names(err) <- ARMS
  sel <- lapply(c("FULL", "B0"), function(a) character(R_REP)); names(sel) <- c("FULL", "B0")

  for (r in seq_len(R_REP)) {
    v1 <- rnorm(N); v2 <- rnorm(N)
    y  <- kg2_predict(THETA_TRUE, v1, v2) + SIGMA * gen(N)

    e_full <- inner_estimate(v1, v2, y, ctrl_full(seed * 1000L + r))
    e_b0   <- inner_estimate(v1, v2, y, ctrl_b0)
    e_pmm2 <- inner_estimate(v1, v2, y, ctrl_pmm2)
    e_pmm3 <- inner_estimate(v1, v2, y, ctrl_pmm3)
    e_lse  <- inner_estimate(v1, v2, y, ctrl_lse)

    err$FULL[r] <- sum((e_full$theta - THETA_TRUE)^2)
    err$B0[r]   <- sum((e_b0$theta   - THETA_TRUE)^2)
    err$PMM2[r] <- sum((e_pmm2$theta - THETA_TRUE)^2)
    err$PMM3[r] <- sum((e_pmm3$theta - THETA_TRUE)^2)
    err$LSE[r]  <- sum((e_lse$theta  - THETA_TRUE)^2)

    sel$FULL[r] <- e_full$method   # realized method (after any PMM->LSE fallback)
    sel$B0[r]   <- e_b0$method
  }

  tgt <- pop_target(gen)
  mse_lse <- mean(err$LSE)
  are <- function(a) round(mse_lse / mean(err[[a]]), 3)

  # PMM activation share = fraction of reps a dispatch arm chose PMM2 or PMM3
  pmm_share <- function(s) round(mean(s %in% c("PMM2", "PMM3")), 2)
  modal     <- function(s) names(sort(table(s), decreasing = TRUE))[1]
  share_str <- function(s) {
    tb <- table(factor(s, levels = c("LSE", "PMM2", "PMM3")))
    fr <- round(as.numeric(tb) / length(s), 2)
    paste(sprintf("%s=%.2f", names(tb), fr), collapse = " ")
  }

  data.frame(
    regime       = name,
    gamma3       = round(tgt$g3, 3),
    gamma4       = round(tgt$g4, 3),
    favored      = tgt$m,
    # selected method (modal) for the two dispatch arms
    sel_FULL     = modal(sel$FULL),
    sel_B0       = modal(sel$B0),
    pmm_share_FULL = pmm_share(sel$FULL),
    pmm_share_B0   = pmm_share(sel$B0),
    # ARE = MSE(LSE)/MSE(arm) for every arm
    ARE_FULL     = are("FULL"),
    ARE_B0       = are("B0"),
    ARE_PMM2     = are("PMM2"),
    ARE_PMM3     = are("PMM3"),
    ARE_LSE      = are("LSE"),     # == 1.000 by construction (sanity)
    detail_FULL  = share_str(sel$FULL),
    detail_B0    = share_str(sel$B0),
    stringsAsFactors = FALSE
  )
}

cat(sprintf("=== E3 ablation (R=%d, n=%d, B=%d) ===\n\n", R_REP, N, B))
t0  <- Sys.time()
res <- do.call(rbind, Map(run_cell, names(gens), gens, seq_along(gens) + 2000L))

# ---- console report --------------------------------------------------------
cat("--- Selected method + PMM activation share (dispatch arms) ---\n")
print(res[, c("regime", "favored", "sel_FULL", "sel_B0",
              "pmm_share_FULL", "pmm_share_B0")], row.names = FALSE)

cat("\n--- Empirical ARE = MSE(LSE)/MSE(arm)  (>1 beats LSE) ---\n")
print(res[, c("regime", "favored", "ARE_FULL", "ARE_B0",
              "ARE_PMM2", "ARE_PMM3", "ARE_LSE")], row.names = FALSE)

cat("\n--- Dispatch detail (per-method selection fraction) ---\n")
print(res[, c("regime", "favored", "detail_FULL", "detail_B0")], row.names = FALSE)

# ---- claim (a): bootstrap suppresses false PMM on heavy-tailed decoys ------
decoys <- res[res$regime %in% c("laplace", "contam"), ]
cat("\n=== CLAIM (a): bootstrap suppresses FALSE PMM activation on decoys ===\n")
for (i in seq_len(nrow(decoys))) {
  cat(sprintf("  %-8s : PMM share  FULL(B=%d)=%.2f  vs  B0=%.2f\n",
              decoys$regime[i], B, decoys$pmm_share_FULL[i], decoys$pmm_share_B0[i]))
}
claim_a <- all(decoys$pmm_share_FULL <= decoys$pmm_share_B0)
cat(sprintf("  => claim (a) holds (FULL <= B0 false-activation on both decoys): %s\n", claim_a))

# ---- claim (b): always-forcing PMM hurts on Gaussian + decoys -------------
gd <- res[res$regime %in% c("Gaussian", "laplace", "contam"), ]
cat("\n=== CLAIM (b): forcing PMM HURTS on Gaussian + decoys vs dispatch ===\n")
for (i in seq_len(nrow(gd))) {
  worst_forced <- min(gd$ARE_PMM2[i], gd$ARE_PMM3[i])
  cat(sprintf("  %-8s : ARE_FULL=%.3f | ARE_PMM2=%.3f ARE_PMM3=%.3f (forced<1 => worse than LSE)\n",
              gd$regime[i], gd$ARE_FULL[i], gd$ARE_PMM2[i], gd$ARE_PMM3[i]))
}
claim_b <- all(pmin(gd$ARE_PMM2, gd$ARE_PMM3) < 1) &&
           all(gd$ARE_FULL >= pmin(gd$ARE_PMM2, gd$ARE_PMM3))
cat(sprintf("  => claim (b) holds (some forced-PMM < LSE, and FULL >= worst forced on each): %s\n", claim_b))

# ---- persist ---------------------------------------------------------------
outdir <- "experiments/revision"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(res, file.path(outdir, "e3_ablation.csv"), row.names = FALSE)
cat(sprintf("\nSaved %s | elapsed %.1fs\n",
            file.path(outdir, "e3_ablation.csv"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
