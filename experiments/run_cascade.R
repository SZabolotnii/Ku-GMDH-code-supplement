#!/usr/bin/env Rscript
# GMDH-level cascade experiment: function recovery and H4 cascade stability
# (experiment-protocol.md 3.2, 3.4, 3.6 H4).
#
# A full GMDH tournament (not a single partial model) is grown to a fixed depth
# on data from a known nonlinear target. We compare the automatic PMM dispatch
# against forced-LSE, ridge-LSE, Huber, and L1 tournaments under Gaussian vs
# strongly skewed noise, and track, by depth L, the best model's:
#   * recovery RMSE     ||f_hat - f_true|| / ||f_true - mean(f_true)||  (test set)
#   * prediction bias   mean(f_hat - f_true)                            (test set)
#   * prediction SD across Monte Carlo repetitions (cascade variance)
#
# Note (honest framing): OLS is unbiased for any mean-zero noise, so any cascade
# "bias" arises from inner-estimator *variance* propagating through the squaring
# compositions (Jensen) and from selection error. The robust signal H4 can show
# is that PMM keeps a lower-variance, lower-RMSE recovery *sustained through
# depth*; we report whatever the numbers say.
#
# Run from the package directory:  Rscript experiments/run_cascade.R [R] [n] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args  <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 40L
N     <- if (length(args) >= 2) as.integer(args[2]) else 800L
B     <- if (length(args) >= 3) as.integer(args[3]) else 150L
N_TEST <- 2000L
L_MAX  <- 4L
D      <- 5L                      # 4 active inputs + 1 decoy (x5)

# Known nonlinear target (degree 2 on x1..x4; x5 is a decoy GMDH must ignore).
f_true <- function(X) {
  1 + 1.5 * X[, 1] - X[, 2] + 0.8 * X[, 1] * X[, 2] +
    0.6 * X[, 1]^2 - 0.4 * X[, 2]^2 + 0.5 * X[, 3] * X[, 4]
}

std <- function(x) (x - mean(x)) / stats::sd(x)
noise_gens <- list(
  Gaussian = function(n) rnorm(n),
  Skew     = function(n) rexp(n) - 1            # gamma3 = 2, gamma4 = 6 -> PMM2
)

ctrl_pmm <- function() gmdh_pmm_control(B = B, L_max = L_MAX, F = 8L, epsilon = -Inf)
ctrl_lse <- function() gmdh_pmm_control(B = 0, L_max = L_MAX, F = 8L, epsilon = -Inf,
                                        force_method = "LSE")
ctrl_ridge_lse <- function() gmdh_pmm_control(B = 0, L_max = L_MAX, F = 8L,
                                              epsilon = -Inf,
                                              force_method = "ridge-LSE",
                                              ridge_lambda = 1e-8)
ctrl_huber <- function() gmdh_pmm_control(B = 0, L_max = L_MAX, F = 8L,
                                          epsilon = -Inf,
                                          force_method = "Huber",
                                          ridge_lambda = 1e-8)
ctrl_l1 <- function() gmdh_pmm_control(B = 0, L_max = L_MAX, F = 8L,
                                       epsilon = -Inf,
                                       force_method = "L1",
                                       ridge_lambda = 1e-8)

# one Monte Carlo replication -> per-depth recovery for both arms
one_rep <- function(gen, seed) {
  set.seed(seed)
  Xtr <- matrix(rnorm(N * D), N, D); ftr <- f_true(Xtr)
  ytr <- ftr + std(gen(N)) * 1.0
  Xte <- matrix(rnorm(N_TEST * D), N_TEST, D); fte <- f_true(Xte)

  fits <- list(PMM = gmdh_pmm(Xtr, ytr, ctrl_pmm()),
               LSE = gmdh_pmm(Xtr, ytr, ctrl_lse()),
               `ridge-LSE` = gmdh_pmm(Xtr, ytr, ctrl_ridge_lse()),
               Huber = gmdh_pmm(Xtr, ytr, ctrl_huber()),
               L1 = gmdh_pmm(Xtr, ytr, ctrl_l1()))
  denom <- sqrt(mean((fte - mean(fte))^2))

  out <- list()
  for (arm in names(fits)) {
    fit <- fits[[arm]]
    for (L in seq_len(min(L_MAX, length(fit$layers)))) {
      yh <- predict(fit, Xte, node = fit$layers[[L]]$best_id)
      out[[length(out) + 1]] <- data.frame(
        arm = arm, depth = L,
        rmse = sqrt(mean((yh - fte)^2)) / denom,
        bias = mean(yh - fte),
        pmm2 = as.integer(fit$layers[[L]]$methods["PMM2"]),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

summarize <- function(label, gen, base_seed) {
  reps <- lapply(seq_len(R_REP), function(r) {
    pr <- try(one_rep(gen, base_seed + r), silent = TRUE)
    if (inherits(pr, "try-error")) NULL else { pr$rep <- r; pr }
  })
  df <- do.call(rbind, Filter(Negate(is.null), reps))
  agg <- do.call(rbind, lapply(split(df, list(df$arm, df$depth)), function(g) {
    data.frame(noise = label, arm = g$arm[1], depth = g$depth[1],
               rmse_mean = round(mean(g$rmse), 3),
               rmse_sd   = round(stats::sd(g$rmse), 3),
               bias_mean = round(mean(g$bias), 4),
               pmm2_frac = round(mean(g$pmm2 > 0, na.rm = TRUE), 2),
               stringsAsFactors = FALSE)
  }))
  agg[order(agg$arm, agg$depth), ]
}

cat(sprintf("=== GMDH-PMM cascade / recovery (R=%d, n=%d, test=%d, L=%d, B=%d) ===\n\n",
            R_REP, N, N_TEST, L_MAX, B))
t0 <- Sys.time()
res <- rbind(summarize("Gaussian", noise_gens$Gaussian, 5000L),
             summarize("Skew",     noise_gens$Skew,     7000L))
print(res, row.names = FALSE)

# --- verdicts --------------------------------------------------------------
deepest <- max(res$depth)
skew_deep <- res[res$noise == "Skew" & res$depth == deepest, ]
r_pmm <- skew_deep$rmse_mean[skew_deep$arm == "PMM"]
r_lse <- skew_deep$rmse_mean[skew_deep$arm == "LSE"]
r_ridge <- skew_deep$rmse_mean[skew_deep$arm == "ridge-LSE"]
r_huber <- skew_deep$rmse_mean[skew_deep$arm == "Huber"]
r_l1 <- skew_deep$rmse_mean[skew_deep$arm == "L1"]
sd_pmm <- skew_deep$rmse_sd[skew_deep$arm == "PMM"]
sd_lse <- skew_deep$rmse_sd[skew_deep$arm == "LSE"]
sd_ridge <- skew_deep$rmse_sd[skew_deep$arm == "ridge-LSE"]
sd_huber <- skew_deep$rmse_sd[skew_deep$arm == "Huber"]
sd_l1 <- skew_deep$rmse_sd[skew_deep$arm == "L1"]

cat(sprintf(
  "\nAt depth %d under Skew:  recovery RMSE  PMM=%.3f vs LSE=%.3f vs ridge-LSE=%.3f vs Huber=%.3f vs L1=%.3f\n",
  deepest, r_pmm, r_lse, r_ridge, r_huber, r_l1))
cat(sprintf(
  "                        RMSE SD       PMM=%.3f vs LSE=%.3f vs ridge-LSE=%.3f vs Huber=%.3f vs L1=%.3f\n",
  sd_pmm, sd_lse, sd_ridge, sd_huber, sd_l1))
cat(sprintf(">>> PMM beats unregularized LSE at depth under skew: %s\n", r_pmm < r_lse))
cat(sprintf(">>> PMM beats ridge-LSE at depth under skew: %s\n", r_pmm < r_ridge))
cat(sprintf(">>> PMM beats Huber at depth under skew: %s\n", r_pmm < r_huber))
cat(sprintf(">>> PMM beats L1 at depth under skew: %s\n", r_pmm < r_l1))

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(res, file.path(outdir, "cascade_recovery.csv"), row.names = FALSE)

# recovery-vs-depth plot (PMM vs LSE, by noise)
png(file.path(outdir, "cascade_recovery.png"), width = 820, height = 460, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (nz in c("Gaussian", "Skew")) {
  sub <- res[res$noise == nz, ]
  yl <- range(sub$rmse_mean)
  plot(NA, xlim = c(1, deepest), ylim = yl, xlab = "GMDH depth", ylab = "recovery RMSE",
       main = nz)
  for (a in c("PMM", "LSE", "ridge-LSE", "Huber", "L1")) {
    s <- sub[sub$arm == a, ]
    lines(s$depth, s$rmse_mean, type = "b", pch = 19,
          col = switch(a, PMM = "firebrick", LSE = "grey45", `ridge-LSE` = "steelblue",
                       Huber = "darkgreen", L1 = "purple4"))
  }
  legend("topright", c("PMM", "LSE", "ridge-LSE", "Huber", "L1"),
         col = c("firebrick", "grey45", "steelblue", "darkgreen", "purple4"),
         pch = 19, lty = 1, bty = "n", cex = 0.8)
}
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s (+ cascade_recovery.png) | elapsed %.1fs\n",
            file.path(outdir, "cascade_recovery.csv"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
