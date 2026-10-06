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

# Every arm sees the SAME train/validation split and the same bootstrap seed.
# (Before 2026-10-06 each arm drew its own split from the global RNG, so arm
# differences mixed estimator effects with split noise; under Gaussian noise
# that alone produced a spurious PMM-over-LSE gap at depth 4.)
ctrl_arm <- function(arm, train, val, seed) {
  base <- list(L_max = L_MAX, F = 8L, epsilon = -Inf,
               train_index = train, val_index = val, boot_seed = seed)
  extra <- switch(arm,
                  PMM = list(B = B),
                  LSE = list(B = 0, force_method = "LSE"),
                  list(B = 0, force_method = arm, ridge_lambda = 1e-8))
  do.call(gmdh_pmm_control, c(base, extra))
}
ARMS <- c("PMM", "LSE", "ridge-LSE", "Huber", "L1")

# one Monte Carlo replication -> per-depth recovery for every arm
one_rep <- function(gen, seed) {
  set.seed(seed)
  Xtr <- matrix(rnorm(N * D), N, D); ftr <- f_true(Xtr)
  ytr <- ftr + std(gen(N)) * 1.0
  Xte <- matrix(rnorm(N_TEST * D), N_TEST, D); fte <- f_true(Xte)
  train <- sort(sample.int(N, floor(0.6 * N)))
  val <- setdiff(seq_len(N), train)

  fits <- stats::setNames(lapply(ARMS, function(a)
    gmdh_pmm(Xtr, ytr, ctrl_arm(a, train, val, seed))), ARMS)
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
        pmm2 = {m <- fit$layers[[L]]$methods; sum(m[names(m) == "PMM2"]) / sum(m)},
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
  df$noise <- label
  raw_all[[label]] <<- df
  agg <- do.call(rbind, lapply(split(df, list(df$arm, df$depth)), function(g) {
    lse <- df[df$arm == "LSE" & df$depth == g$depth[1], ]
    m <- merge(g[, c("rep", "rmse")], lse[, c("rep", "rmse")], by = "rep")
    wp <- if (g$arm[1] == "LSE") NA_real_ else
      suppressWarnings(stats::wilcox.test(m$rmse.x, m$rmse.y, paired = TRUE)$p.value)
    data.frame(noise = label, arm = g$arm[1], depth = g$depth[1], reps = nrow(g),
               rmse_median = round(stats::median(g$rmse), 3),
               rmse_q25  = round(stats::quantile(g$rmse, 0.25, names = FALSE), 3),
               rmse_q75  = round(stats::quantile(g$rmse, 0.75, names = FALSE), 3),
               rmse_mean = round(mean(g$rmse), 3),
               rmse_sd   = round(stats::sd(g$rmse), 3),
               rmse_max  = round(max(g$rmse), 3),
               n_blowup  = sum(g$rmse > 0.5),
               ratio_vs_lse_median = round(stats::median(m$rmse.x / m$rmse.y), 3),
               wins_vs_lse = sum(m$rmse.x < m$rmse.y),
               wilcox_p_vs_lse = signif(wp, 3),
               bias_mean = round(mean(g$bias), 4),
               pmm2_share = round(mean(g$pmm2, na.rm = TRUE), 2),
               stringsAsFactors = FALSE)
  }))
  agg[order(agg$arm, agg$depth), ]
}

cat(sprintf("=== GMDH-PMM cascade / recovery (R=%d, n=%d, test=%d, L=%d, B=%d) ===\n\n",
            R_REP, N, N_TEST, L_MAX, B))
t0 <- Sys.time()
raw_all <- list()
res <- rbind(summarize("Gaussian", noise_gens$Gaussian, 5000L),
             summarize("Skew",     noise_gens$Skew,     7000L))
print(res, row.names = FALSE)

# --- verdicts (medians + paired tests; depth-4 means are driven by blow-ups) --
for (nz in unique(res$noise)) for (L in sort(unique(res$depth))) {
  sub <- res[res$noise == nz & res$depth == L & res$arm != "LSE", ]
  for (i in seq_len(nrow(sub)))
    cat(sprintf("%-8s depth %d  %-9s median ratio vs LSE %.3f  wins %d/%d  Wilcoxon p=%.3g\n",
                nz, L, sub$arm[i], sub$ratio_vs_lse_median[i], sub$wins_vs_lse[i],
                sub$reps[i], sub$wilcox_p_vs_lse[i]))
}

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
# _v2: common-split protocol (2026-10-06). cascade_recovery.csv is the pre-fix
# artifact behind Kybernetika 7414 Tab. cascade and is kept for provenance.
utils::write.csv(res, file.path(outdir, "cascade_recovery_v2.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, raw_all), file.path(outdir, "cascade_recovery_v2_raw.csv"),
                 row.names = FALSE)

# recovery-vs-depth plot (PMM vs LSE, by noise)
deepest <- max(res$depth)
png(file.path(outdir, "cascade_recovery_v2.png"), width = 820, height = 460, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (nz in c("Gaussian", "Skew")) {
  sub <- res[res$noise == nz, ]
  yl <- range(sub$rmse_median)
  plot(NA, xlim = c(1, deepest), ylim = yl, xlab = "GMDH depth", ylab = "median recovery NRMSE",
       main = nz)
  for (a in c("PMM", "LSE", "ridge-LSE", "Huber", "L1")) {
    s <- sub[sub$arm == a, ]
    lines(s$depth, s$rmse_median, type = "b", pch = 19,
          col = switch(a, PMM = "firebrick", LSE = "grey45", `ridge-LSE` = "steelblue",
                       Huber = "darkgreen", L1 = "purple4"))
  }
  legend("topright", c("PMM", "LSE", "ridge-LSE", "Huber", "L1"),
         col = c("firebrick", "grey45", "steelblue", "darkgreen", "purple4"),
         pch = 19, lty = 1, bty = "n", cex = 0.8)
}
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s (+ _raw.csv, .png) | elapsed %.1fs\n",
            file.path(outdir, "cascade_recovery_v2.csv"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
