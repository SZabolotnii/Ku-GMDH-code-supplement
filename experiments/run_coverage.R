#!/usr/bin/env Rscript
# H1 coverage experiment for Paper 1.
#
# This is a synthetic inferential check before the real-world C5 stage. It asks:
# do PMM-based coefficient intervals preserve nominal 95% coefficient coverage
# while being narrower than LSE/robust baselines in PMM-favorable regimes?
#
# We deliberately report both coverage and width. A method that narrows
# intervals by under-covering is not a success.
#
# Run from the package directory:
#   Rscript experiments/run_coverage.R [R] [n] [B_ci] [B_diag]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
R_REP  <- if (length(args) >= 1) as.integer(args[1]) else 60L
N      <- if (length(args) >= 2) as.integer(args[2]) else 500L
B_CI   <- if (length(args) >= 3) as.integer(args[3]) else 80L
B_DIAG <- if (length(args) >= 4) as.integer(args[4]) else 80L
LEVEL  <- 0.95
THETA_TRUE <- c("(Intercept)" = 0.5, b1 = 1.5, b2 = -1.0, b12 = 0.8, b11 = 0.6, b22 = -0.4)
COEF_ORDER <- names(THETA_TRUE)

std <- function(x) (x - mean(x)) / stats::sd(x)
gens <- list(
  exp = function(n) rexp(n) - 1,
  uniform = function(n) runif(n, -sqrt(3), sqrt(3))
)

ctrl_for <- function(method) {
  switch(method,
         LSE = gmdh_pmm_control(B = 0, force_method = "LSE"),
         PMM = gmdh_pmm_control(B = B_DIAG, force_method = "auto"),
         Huber = gmdh_pmm_control(B = 0, force_method = "Huber"),
         L1 = gmdh_pmm_control(B = 0, force_method = "L1"),
         stop("Unknown method: ", method))
}

fit_theta <- function(v1, v2, y, method) {
  inner_estimate(v1, v2, y, ctrl_for(method))$theta[COEF_ORDER]
}

fit_lse_normal_ci <- function(v1, v2, y) {
  d <- kg2_design(v1, v2, y)
  fit <- stats::lm(y ~ b1 + b2 + b12 + b11 + b22, data = d)
  cf <- stats::coef(fit)[COEF_ORDER]
  cf[is.na(cf)] <- 0
  sm <- summary(fit)$coefficients
  se <- sm[COEF_ORDER, "Std. Error"]
  se[!is.finite(se)] <- Inf
  crit <- stats::qt((1 + LEVEL) / 2, stats::df.residual(fit))
  list(theta = cf, lower = cf - crit * se, upper = cf + crit * se)
}

bootstrap_ci <- function(v1, v2, y, method, seed) {
  set.seed(seed)
  n <- length(y)
  theta_hat <- fit_theta(v1, v2, y, method)
  boot <- matrix(NA_real_, nrow = B_CI, ncol = length(theta_hat),
                 dimnames = list(NULL, names(theta_hat)))
  for (b in seq_len(B_CI)) {
    idx <- sample.int(n, n, replace = TRUE)
    th <- try(fit_theta(v1[idx], v2[idx], y[idx], method), silent = TRUE)
    if (!inherits(th, "try-error") && all(is.finite(th))) boot[b, ] <- th
  }
  boot_sd <- apply(boot, 2, stats::sd, na.rm = TRUE)
  z <- stats::qnorm((1 + LEVEL) / 2)
  lower_norm <- theta_hat - z * boot_sd
  upper_norm <- theta_hat + z * boot_sd
  qs <- apply(boot, 2, stats::quantile,
              probs = c((1 - LEVEL) / 2, (1 + LEVEL) / 2),
              na.rm = TRUE, names = FALSE)
  list(theta = theta_hat,
       normal = list(lower = qsafe(lower_norm), upper = qsafe(upper_norm)),
       percentile = list(lower = qsafe(qs[1, ]), upper = qsafe(qs[2, ])),
       boot_ok = mean(stats::complete.cases(boot)))
}

qsafe <- function(x) {
  x[!is.finite(x)] <- NA_real_
  x
}

coverage_row <- function(scenario, rep_id, method, ci_type, theta, lower, upper, boot_ok = NA_real_) {
  covered <- lower <= THETA_TRUE & THETA_TRUE <= upper
  width <- upper - lower
  data.frame(
    scenario = scenario,
    rep = rep_id,
    method = method,
    ci_type = ci_type,
    coef_coverage = mean(covered, na.rm = TRUE),
    joint_coverage = as.numeric(all(covered, na.rm = TRUE)),
    mean_width = mean(width, na.rm = TRUE),
    theta_mse = mean((theta - THETA_TRUE)^2),
    boot_ok = boot_ok,
    stringsAsFactors = FALSE
  )
}

one_rep <- function(scenario, gen, rep_id, seed) {
  set.seed(seed)
  v1 <- rnorm(N); v2 <- rnorm(N)
  y <- kg2_predict(THETA_TRUE, v1, v2) + gen(N)

  out <- list()
  lse_norm <- fit_lse_normal_ci(v1, v2, y)
  out[[length(out) + 1]] <- coverage_row(scenario, rep_id, "LSE", "classical-normal",
                                         lse_norm$theta, lse_norm$lower, lse_norm$upper)

  for (method in c("LSE", "PMM", "Huber", "L1")) {
    ci <- bootstrap_ci(v1, v2, y, method, seed + match(method, c("LSE", "PMM", "Huber", "L1")) * 100000L)
    out[[length(out) + 1]] <- coverage_row(scenario, rep_id, method, "bootstrap-normal",
                                           ci$theta, ci$normal$lower, ci$normal$upper, ci$boot_ok)
    out[[length(out) + 1]] <- coverage_row(scenario, rep_id, method, "bootstrap-percentile",
                                           ci$theta, ci$percentile$lower, ci$percentile$upper, ci$boot_ok)
  }
  do.call(rbind, out)
}

summarize <- function(raw) {
  parts <- split(raw, list(raw$scenario, raw$method, raw$ci_type), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(g) {
    data.frame(
      scenario = g$scenario[1],
      method = g$method[1],
      ci_type = g$ci_type[1],
      coef_coverage = round(mean(g$coef_coverage, na.rm = TRUE), 3),
      joint_coverage = round(mean(g$joint_coverage, na.rm = TRUE), 3),
      mean_width = round(mean(g$mean_width, na.rm = TRUE), 3),
      theta_rmse = round(sqrt(mean(g$theta_mse, na.rm = TRUE)), 3),
      boot_ok = round(mean(g$boot_ok, na.rm = TRUE), 3),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$scenario, out$method, out$ci_type), ]
}

cat(sprintf("=== GMDH-PMM H1 coverage benchmark (R=%d, n=%d, B_ci=%d, B_diag=%d) ===\n\n",
            R_REP, N, B_CI, B_DIAG))
t0 <- Sys.time()
raw <- do.call(rbind, unlist(Map(function(name, gen, offset) {
  lapply(seq_len(R_REP), function(r) one_rep(name, gen, r, 9000L + offset * 1000L + r))
}, names(gens), gens, seq_along(gens)), recursive = FALSE))
res <- summarize(raw)
print(res, row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(raw, file.path(outdir, "coverage_h1_raw.csv"), row.names = FALSE)
utils::write.csv(res, file.path(outdir, "coverage_h1.csv"), row.names = FALSE)

png(file.path(outdir, "coverage_h1.png"), width = 920, height = 520, res = 110)
op <- par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
for (sc in names(gens)) {
  sub <- res[res$scenario == sc & res$ci_type != "bootstrap-percentile", ]
  lab <- paste(sub$method, sub$ci_type, sep = "\n")
  cols <- ifelse(sub$method == "PMM", "firebrick",
                 ifelse(sub$method == "LSE", "grey45",
                        ifelse(sub$method == "Huber", "darkgreen", "purple4")))
  bp <- barplot(sub$coef_coverage, names.arg = lab, ylim = c(0, 1.05),
                las = 2, col = cols, ylab = "coefficient-wise coverage",
                main = sc, cex.names = 0.65)
  abline(h = LEVEL, lty = 2, col = "grey20")
  text(bp, pmin(sub$coef_coverage + 0.04, 1.03), sub$mean_width, cex = 0.65)
  mtext("numbers above bars = mean interval width", side = 3, line = 0.2, cex = 0.65)
}
par(op); invisible(dev.off())

cat(sprintf("\nSaved coverage_h1.csv, coverage_h1_raw.csv, coverage_h1.png | elapsed %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
