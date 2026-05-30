#!/usr/bin/env Rscript
# EXPERIMENT E6 -- MINIMAL-MOMENT GMM BASELINE (reviewer Q6).
#
# Question: is PMM2 just a re-branding of GMM on two moments, or a competitor?
# Theory (Kunchenko / Godambe optimal estimating functions): PMM2 IS the
# optimally-weighted estimator using exactly the score from the mean moment and
# the squared-residual moment. So a *two-step efficient* GMM using the moment
# conditions
#       g(theta) = (1/n) [ Z^T e ; Z^T e^2 ],   e = y - Z theta
# should reproduce PMM2 closely -- confirming PMM is the *economical* member of
# the GMM family, not a rival to it.
#
# Design: canonical KG-2 linear-in-parameters regression (same as
# experiments/run_synthetic.R), y = Z theta + e, with
#   Z = [1, v1, v2, v1*v2, v1^2, v2^2]   (n x 6, theta length 6).
# Regimes: exp-skew (rexp-1) and lognorm-skew (std lognormal) -- both PMM2-favorable.
# Metrics: coefficient MSE = mean over reps of ||theta_hat - theta_true||^2,
#          and test NRMSE on an independent test fold.
# Estimators: GMM(e, e^2) [two-step efficient], PMM2 (force_method="PMM2"), LSE.
#
# Run from the package dir:
#   Rscript experiments/revision/e6_gmm.R [REPS] [N]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args  <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 50L
N     <- if (length(args) >= 2) as.integer(args[2]) else 2000L
SIGMA <- 1.0
THETA_TRUE <- c(0.5, 1.5, -1.0, 0.8, 0.6, -0.4)  # canonical KG-2 order

std <- function(x) (x - mean(x)) / stats::sd(x)
gens <- list(
  exp     = function(n) rexp(n) - 1,            # strong + skew, gamma3 ~ 2
  lognorm = function(n) std(rlnorm(n, 0, 0.5))  # A3-like skew, gamma3 ~ 1.75
)

# Canonical KG-2 design matrix (n x 6), columns in THETA_TRUE order.
Zmat <- function(v1, v2) cbind(1, v1, v2, v1 * v2, v1^2, v2^2)

# ---------------------------------------------------------------------------
# Two-moment two-step efficient GMM for y = Z theta + e.
#
# CRUCIAL: the squared-residual moment must be CENTRED. The raw moment
# E[z_i e_i^2] equals z_i * sigma^2 at the truth, which is NOT zero for any
# column of Z with non-zero mean (e.g. the intercept). Stacking the raw e^2
# moment is therefore MISSPECIFIED -- it forces mean(e^2) -> 0 and wrecks the
# fit. PMM2's optimal estimating function uses the centred second moment, so
# the faithful "two-moment GMM" introduces the residual variance s2 = E[e^2] as
# a nuisance parameter and uses the moment block z_i (e_i^2 - s2).
#
# Parameters: beta = (theta, s2), length p+1.
# Per-observation moment vector (length 2p):
#       g_i(beta) = [ z_i e_i ; z_i (e_i^2 - s2) ],   e_i = y_i - z_i' theta.
# (The intercept column of the second block, z=1, supplies mean(e^2)=s2, so s2
#  is identified; the system is just-identified: 2p moments, p+1 params, and the
#  extra p-1 over-identifying conditions are exactly the skew-informative ones
#  that make the optimal weighting outperform LSE -- the same content as PMM2.)
#
# Sample moment gbar(beta) = (1/n) sum_i g_i(beta).  Objective Q = gbar' W gbar.
# Step 1: W = I_{2p}.   Step 2: W = Omega^{-1}, Omega = (1/n) sum g_i g_i' at
# the step-1 estimate (the optimal GMM weighting = inverse moment covariance).
#
# Analytic Jacobian (de_i/dtheta = -z_i, d/ds2 only hits the 2nd block):
#   d/dtheta [ z_i e_i ]        = -z_i z_i'        ;  d/ds2 = 0
#   d/dtheta [ z_i (e_i^2-s2) ] = -2 e_i z_i z_i'  ;  d/ds2 = -z_i
# so  G(beta) = (1/n) [ -Z'Z ,        0      ;
#                       -2 Z'diag(e)Z , -Z'1 ]    (2p x (p+1)).
# Gauss-Newton with step-halving line search on the GMM objective.
# ---------------------------------------------------------------------------
gmm_two_moment <- function(Z, y, theta0, s2_0 = NULL, max_iter = 200L, tol = 1e-10) {
  p  <- ncol(Z)
  n  <- nrow(Z)
  ZtZ <- crossprod(Z)
  Zt1 <- colSums(Z)                         # Z' 1  (length p)
  if (is.null(s2_0)) s2_0 <- mean((y - Z %*% theta0)^2)
  beta0 <- c(theta0, s2_0)

  moments <- function(beta) {
    theta <- beta[1:p]; s2 <- beta[p + 1]
    e  <- as.numeric(y - Z %*% theta)
    gi <- cbind(Z * e, Z * (e^2 - s2))      # n x 2p per-obs moments
    list(e = e, gi = gi, gbar = colMeans(gi))
  }
  jacobian <- function(e) {
    top <- cbind(-ZtZ / n, rep(0, p))                       # p x (p+1)
    bot <- cbind(-2 * crossprod(Z, Z * e) / n, -Zt1 / n)    # p x (p+1)
    rbind(top, bot)                                         # 2p x (p+1)
  }
  solveQ <- function(beta, W) {            # Gauss-Newton minimization of g'Wg
    b <- beta
    for (it in seq_len(max_iter)) {
      m  <- moments(b)
      G  <- jacobian(m$e)
      Q0 <- as.numeric(t(m$gbar) %*% W %*% m$gbar)
      GtW <- crossprod(G, W)               # (p+1) x 2p
      A   <- GtW %*% G                      # (p+1) x (p+1)
      rhs <- GtW %*% m$gbar                 # (p+1) x 1
      step <- tryCatch(as.numeric(solve(A, rhs)),
                       error = function(z) as.numeric(MASS::ginv(A) %*% rhs))
      # step-halving line search on the GMM objective
      lam <- 1; improved <- FALSE
      for (h in 0:40) {
        b_try <- b - lam * step
        m_try <- moments(b_try)
        Q1 <- as.numeric(t(m_try$gbar) %*% W %*% m_try$gbar)
        if (is.finite(Q1) && Q1 < Q0) { b <- b_try; improved <- TRUE; break }
        lam <- lam / 2
      }
      if (!improved) break
      if (sqrt(sum((lam * step)^2)) < tol * (1 + sqrt(sum(b^2)))) break
    }
    b
  }

  # Step 1: identity weighting.
  W1    <- diag(2 * p)
  beta1 <- solveQ(beta0, W1)

  # Step 2: efficient weighting from the step-1 moment covariance.
  m1    <- moments(beta1)
  Omega <- crossprod(m1$gi) / n                       # 2p x 2p (E[g g'] at beta1)
  Omega <- Omega + diag(1e-8, 2 * p)                  # ridge for stability
  W2    <- tryCatch(solve(Omega), error = function(z) MASS::ginv(Omega))
  beta_eff <- solveQ(beta1, W2)

  list(theta = beta_eff[1:p], s2 = beta_eff[p + 1],
       theta_step1 = beta1[1:p])
}

ctrl_lse  <- gmdh_pmm_control(B = 0, force_method = "LSE")
ctrl_pmm2 <- gmdh_pmm_control(B = 0, force_method = "PMM2")

nrmse <- function(pred, y) sqrt(mean((pred - y)^2)) / stats::sd(y)

run_regime <- function(name, gen, seed0) {
  err_lse <- err_pmm2 <- err_gmm <- numeric(REPS)
  nr_lse  <- nr_pmm2  <- nr_gmm  <- numeric(REPS)
  d_gmm_pmm2 <- numeric(REPS)
  conv_gmm <- 0L
  for (r in seq_len(REPS)) {
    set.seed(seed0 + r)
    v1 <- rnorm(N); v2 <- rnorm(N)
    y  <- kg2_predict(THETA_TRUE, v1, v2) + SIGMA * gen(N)
    v1t <- rnorm(N); v2t <- rnorm(N)
    yt  <- kg2_predict(THETA_TRUE, v1t, v2t) + SIGMA * gen(N)
    Zt  <- Zmat(v1t, v2t)

    th_lse  <- inner_estimate(v1, v2, y, ctrl_lse)$theta
    th_pmm2 <- inner_estimate(v1, v2, y, ctrl_pmm2)$theta
    g       <- gmm_two_moment(Zmat(v1, v2), y, theta0 = as.numeric(th_lse))
    th_gmm  <- g$theta
    if (all(is.finite(th_gmm))) conv_gmm <- conv_gmm + 1L

    err_lse[r]  <- sum((th_lse  - THETA_TRUE)^2)
    err_pmm2[r] <- sum((th_pmm2 - THETA_TRUE)^2)
    err_gmm[r]  <- sum((th_gmm  - THETA_TRUE)^2)
    d_gmm_pmm2[r] <- sqrt(sum((as.numeric(th_gmm) - as.numeric(th_pmm2))^2))

    nr_lse[r]  <- nrmse(as.numeric(Zt %*% th_lse),  yt)
    nr_pmm2[r] <- nrmse(as.numeric(Zt %*% th_pmm2), yt)
    nr_gmm[r]  <- nrmse(as.numeric(Zt %*% th_gmm),  yt)
  }
  data.frame(
    regime        = name,
    n             = N,
    reps          = REPS,
    coefMSE_LSE   = round(mean(err_lse),  5),
    coefMSE_PMM2  = round(mean(err_pmm2), 5),
    coefMSE_GMM   = round(mean(err_gmm),  5),
    NRMSE_LSE     = round(mean(nr_lse),  5),
    NRMSE_PMM2    = round(mean(nr_pmm2), 5),
    NRMSE_GMM     = round(mean(nr_gmm),  5),
    ratio_GMM_PMM2_coefMSE = round(mean(err_gmm) / mean(err_pmm2), 4),
    ratio_GMM_PMM2_NRMSE   = round(mean(nr_gmm)  / mean(nr_pmm2),  4),
    ARE_PMM2_vs_LSE = round(mean(err_lse) / mean(err_pmm2), 3),
    ARE_GMM_vs_LSE  = round(mean(err_lse) / mean(err_gmm),  3),
    mean_coef_dist_GMM_PMM2 = round(mean(d_gmm_pmm2), 5),
    gmm_finite_frac = conv_gmm / REPS,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("=== E6: minimal-moment GMM(e, e^2) vs PMM2 vs LSE  (reps=%d, n=%d) ===\n\n",
            REPS, N))
t0 <- Sys.time()
res <- do.call(rbind, Map(run_regime, names(gens), gens, c(20100L, 20200L)))

# Print in two readable blocks.
cat("--- Coefficient MSE  (mean ||theta_hat - theta_true||^2) ---\n")
print(res[, c("regime", "coefMSE_LSE", "coefMSE_PMM2", "coefMSE_GMM",
              "ratio_GMM_PMM2_coefMSE", "ARE_PMM2_vs_LSE", "ARE_GMM_vs_LSE")],
      row.names = FALSE)
cat("\n--- Test NRMSE ---\n")
print(res[, c("regime", "NRMSE_LSE", "NRMSE_PMM2", "NRMSE_GMM",
              "ratio_GMM_PMM2_NRMSE")], row.names = FALSE)
cat("\n--- GMM<->PMM2 proximity ---\n")
print(res[, c("regime", "mean_coef_dist_GMM_PMM2", "gmm_finite_frac")],
      row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(res, file.path(outdir, "e6_gmm.csv"), row.names = FALSE)
# also drop a copy next to the script
utils::write.csv(res, "experiments/revision/e6_gmm.csv", row.names = FALSE)

cat(sprintf("\nSaved experiments/revision/e6_gmm.csv (+ results/e6_gmm.csv) | elapsed %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
