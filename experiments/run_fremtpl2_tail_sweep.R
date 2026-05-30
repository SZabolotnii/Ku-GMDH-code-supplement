#!/usr/bin/env Rscript
# C5 freMTPL2 TailReserve weight sweep + tail calibration.
#
# This experiment follows the reserve-aware criterion pass. It keeps the
# external criterion in reserve-aware mode and sweeps the top-decile penalty
# weight. Final reporting metrics are unchanged across all arms.
#
# Run from the package directory:
#   Rscript experiments/run_fremtpl2_tail_sweep.R [R] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 8L
B <- if (length(args) >= 2) as.integer(args[2]) else 30L
TEST_FRAC <- 0.25
L_MAX <- 3L
F_KEEP <- 6L
TAIL_WEIGHTS <- c(1, 2, 4, 8, 16)

find_root <- function() {
  here <- normalizePath(getwd(), mustWork = TRUE)
  cur <- here
  repeat {
    if (file.exists(file.path(cur, "ROADMAP.md")) &&
        file.exists(file.path(cur, "paper-1-gmdh-pmm", "code", "DESCRIPTION"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Cannot locate repository root from ", here)
    cur <- parent
  }
}

model_x <- function(formula, data) {
  X <- stats::model.matrix(formula, data = data)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  finite <- apply(X, 2, function(z) all(is.finite(z)))
  varied <- apply(X, 2, function(z) stats::sd(z) > 1e-12)
  X <- X[, finite & varied, drop = FALSE]
  storage.mode(X) <- "double"
  X
}

load_fremtpl2_raw <- function(root) {
  path <- file.path(root, "shared", "datasets", "external", "processed",
                    "fremtpl2_severity_sample.csv")
  if (!file.exists(path)) {
    stop("Missing freMTPL2 prepared sample: ", path,
         "\nRun: Rscript shared/datasets/prepare_fremtpl2.R 4000")
  }
  d <- stats::na.omit(utils::read.csv(path, stringsAsFactors = TRUE))
  form <- ClaimAmount ~ Exposure + VehPower + VehAge + DrivAge + BonusMalus +
    log_density + Area + VehGas
  mf <- stats::model.frame(form, data = d, na.action = stats::na.omit)
  list(
    id = "fremtpl2_severity_raw",
    y = stats::model.response(mf),
    X = model_x(form, mf)
  )
}

ctrl_for <- function(estimator, tail_weight, seed) {
  base <- list(
    L_max = L_MAX, F = F_KEEP, epsilon = -Inf, seed = seed,
    criterion = "reserve-aware",
    reserve_weight = 1, tail_weight = tail_weight,
    overreserve_weight = 0.25, tail_prob = 0.9,
    max_iter = 60L
  )
  switch(estimator,
         `PMM-auto` = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "auto"))),
         `PMM-relaxed` = do.call(gmdh_pmm_control, c(base, list(
           B = B, force_method = "auto", alpha = 0.10,
           skew_min = 0.10, skew_strong = 0.50, g2_threshold = 0.99
         ))),
         LSE = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "LSE"))),
         Huber = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "Huber"))),
         L1 = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "L1"))),
         stop("Unknown estimator: ", estimator))
}

standardize_split <- function(X, train, test) {
  mu <- colMeans(X[train, , drop = FALSE])
  sig <- apply(X[train, , drop = FALSE], 2, stats::sd)
  sig[!is.finite(sig) | sig <= 1e-12] <- 1
  list(
    X_train = sweep(sweep(X[train, , drop = FALSE], 2, mu, "-"), 2, sig, "/"),
    X_test = sweep(sweep(X[test, , drop = FALSE], 2, mu, "-"), 2, sig, "/")
  )
}

method_share <- function(fit, method) {
  tabs <- lapply(fit$layers, function(z) z$methods)
  num <- sum(vapply(tabs, function(t) {
    val <- unname(t[method])
    if (!length(val) || is.na(val)) 0 else val
  }, numeric(1)))
  den <- sum(vapply(tabs, sum, numeric(1)))
  if (den <= 0) NA_real_ else num / den
}

metrics <- function(y, pred) {
  e <- pred - y
  y_scale <- stats::sd(y)
  if (!is.finite(y_scale) || y_scale <= 1e-12) y_scale <- 1
  tail <- y >= as.numeric(stats::quantile(y, 0.90, na.rm = TRUE, names = FALSE))
  data.frame(
    nrmse = sqrt(mean(e^2)) / y_scale,
    mae = mean(abs(e)),
    bias = mean(e),
    reserve_error = sum(pred) / sum(y) - 1,
    abs_reserve_error = abs(sum(pred) / sum(y) - 1),
    under_rate = mean(pred < y),
    tail_mae = mean(abs(e[tail])),
    tail_bias = mean(e[tail]),
    tail_reserve_error = sum(pred[tail]) / sum(y[tail]) - 1,
    abs_tail_reserve_error = abs(sum(pred[tail]) / sum(y[tail]) - 1),
    tail_under_rate = mean(pred[tail] < y[tail]),
    stringsAsFactors = FALSE
  )
}

calibration_bins <- function(y, pred, rep_id, tail_weight, estimator) {
  n <- length(y)
  ord <- order(y)
  bin <- integer(n)
  bin[ord] <- pmin(10L, ceiling(seq_len(n) * 10 / n))
  rows <- lapply(seq_len(10L), function(b) {
    idx <- bin == b
    observed_sum <- sum(y[idx])
    predicted_sum <- sum(pred[idx])
    data.frame(
      rep = rep_id,
      tail_weight = tail_weight,
      estimator = estimator,
      claim_decile = b,
      n = sum(idx),
      observed_sum = observed_sum,
      predicted_sum = predicted_sum,
      reserve_error = if (abs(observed_sum) < 1e-12) 0 else predicted_sum / observed_sum - 1,
      observed_mean = mean(y[idx]),
      predicted_mean = mean(pred[idx]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

one_rep <- function(dat, rep_id, seed) {
  set.seed(seed)
  n <- length(dat$y)
  n_test <- max(10L, ceiling(TEST_FRAC * n))
  test <- sort(sample.int(n, n_test))
  train <- setdiff(seq_len(n), test)
  xs <- standardize_split(dat$X, train, test)
  y_train <- dat$y[train]
  y_test <- dat$y[test]
  estimators <- c("PMM-auto", "PMM-relaxed", "LSE", "Huber", "L1")

  rows <- list()
  cal_rows <- list()
  for (tail_weight in TAIL_WEIGHTS) {
    for (estimator in estimators) {
      fit <- try(gmdh_pmm(xs$X_train, y_train,
                          ctrl_for(estimator, tail_weight, seed)),
                 silent = TRUE)
      if (inherits(fit, "try-error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          rep = rep_id, tail_weight = tail_weight, estimator = estimator, ok = FALSE,
          nrmse = NA_real_, mae = NA_real_, bias = NA_real_,
          reserve_error = NA_real_, abs_reserve_error = NA_real_,
          under_rate = NA_real_, tail_mae = NA_real_, tail_bias = NA_real_,
          tail_reserve_error = NA_real_, abs_tail_reserve_error = NA_real_,
          tail_under_rate = NA_real_, pmm2_share = NA_real_, layers = NA_integer_,
          stringsAsFactors = FALSE)
        next
      }
      pred <- predict(fit, xs$X_test)
      m <- metrics(y_test, pred)
      rows[[length(rows) + 1L]] <- data.frame(
        rep = rep_id, tail_weight = tail_weight, estimator = estimator, ok = TRUE,
        m,
        pmm2_share = method_share(fit, "PMM2"),
        layers = length(fit$layers),
        stringsAsFactors = FALSE
      )
      cal_rows[[length(cal_rows) + 1L]] <- calibration_bins(
        y_test, pred, rep_id, tail_weight, estimator
      )
    }
  }
  raw <- do.call(rbind, rows)
  raw$nrmse_winner <- ave(raw$nrmse, raw$rep, raw$tail_weight,
                          FUN = function(z) z == min(z, na.rm = TRUE)) > 0
  raw$reserve_winner <- ave(raw$abs_reserve_error, raw$rep, raw$tail_weight,
                            FUN = function(z) z == min(z, na.rm = TRUE)) > 0
  raw$tail_reserve_winner <- ave(raw$abs_tail_reserve_error, raw$rep, raw$tail_weight,
                                 FUN = function(z) z == min(z, na.rm = TRUE)) > 0
  list(raw = raw, calibration = do.call(rbind, cal_rows))
}

summarize <- function(raw) {
  ok <- raw[raw$ok, ]
  parts <- split(ok, list(ok$tail_weight, ok$estimator), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(g) {
    data.frame(
      tail_weight = g$tail_weight[1],
      estimator = g$estimator[1],
      nrmse_mean = round(mean(g$nrmse), 3),
      mae_mean = round(mean(g$mae), 1),
      bias_mean = round(mean(g$bias), 1),
      reserve_error_mean = round(mean(g$reserve_error), 3),
      abs_reserve_error_mean = round(mean(g$abs_reserve_error), 3),
      tail_reserve_error_mean = round(mean(g$tail_reserve_error), 3),
      abs_tail_reserve_error_mean = round(mean(g$abs_tail_reserve_error), 3),
      nrmse_win_rate = round(mean(g$nrmse_winner), 3),
      reserve_win_rate = round(mean(g$reserve_winner), 3),
      tail_reserve_win_rate = round(mean(g$tail_reserve_winner), 3),
      pmm2_share = round(mean(g$pmm2_share), 3),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$tail_weight, out$nrmse_mean), ]
}

summarize_calibration <- function(cal) {
  parts <- split(cal, list(cal$tail_weight, cal$estimator, cal$claim_decile), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(g) {
    data.frame(
      tail_weight = g$tail_weight[1],
      estimator = g$estimator[1],
      claim_decile = g$claim_decile[1],
      observed_sum_mean = round(mean(g$observed_sum), 1),
      predicted_sum_mean = round(mean(g$predicted_sum), 1),
      reserve_error_mean = round(mean(g$reserve_error), 3),
      observed_mean = round(mean(g$observed_mean), 1),
      predicted_mean = round(mean(g$predicted_mean), 1),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$tail_weight, out$estimator, out$claim_decile), ]
}

cat(sprintf("=== freMTPL2 TailReserve sweep (R=%d, B=%d) ===\n\n", R_REP, B))
t0 <- Sys.time()
root <- find_root()
dat <- load_fremtpl2_raw(root)
rep_out <- lapply(seq_len(R_REP), function(r) {
  one_rep(dat, r, 104000L + r)
})
raw <- do.call(rbind, lapply(rep_out, `[[`, "raw"))
cal <- do.call(rbind, lapply(rep_out, `[[`, "calibration"))
res <- summarize(raw)
cal_res <- summarize_calibration(cal)

print(res, row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(raw, file.path(outdir, "fremtpl2_tail_sweep_raw.csv"), row.names = FALSE)
utils::write.csv(res, file.path(outdir, "fremtpl2_tail_sweep_summary.csv"), row.names = FALSE)
utils::write.csv(cal, file.path(outdir, "fremtpl2_tail_calibration_raw.csv"), row.names = FALSE)
utils::write.csv(cal_res, file.path(outdir, "fremtpl2_tail_calibration_summary.csv"), row.names = FALSE)

png(file.path(outdir, "fremtpl2_tail_sweep.png"), width = 1120, height = 780, res = 115)
op <- par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
estimators <- c("PMM-auto", "PMM-relaxed", "LSE", "Huber", "L1")
cols <- c("firebrick", "tomato", "grey45", "darkgreen", "seagreen")
plot_metric <- function(metric, ylab, main) {
  plot(NA, xlim = range(TAIL_WEIGHTS), ylim = range(res[[metric]], na.rm = TRUE),
       log = "x", xaxt = "n", xlab = "tail_weight", ylab = ylab, main = main)
  axis(1, at = TAIL_WEIGHTS, labels = TAIL_WEIGHTS)
  for (i in seq_along(estimators)) {
    z <- res[res$estimator == estimators[i], ]
    z <- z[order(z$tail_weight), ]
    lines(z$tail_weight, z[[metric]], type = "b", pch = 16, col = cols[i])
  }
  legend("topright", legend = estimators, col = cols, pch = 16, lty = 1, cex = 0.75)
}
plot_metric("nrmse_mean", "mean", "NRMSE")
plot_metric("abs_reserve_error_mean", "absolute error", "Aggregate reserve")
plot_metric("abs_tail_reserve_error_mean", "absolute error", "Top-decile reserve")
plot_metric("mae_mean", "mean", "MAE")
par(op); invisible(dev.off())

png(file.path(outdir, "fremtpl2_tail_calibration.png"), width = 1120, height = 720, res = 115)
op <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
for (tw in c(4, 16)) {
  z <- cal_res[cal_res$tail_weight == tw, ]
  plot(NA, xlim = c(1, 10), ylim = range(z$reserve_error_mean, na.rm = TRUE),
       xlab = "observed claim decile", ylab = "sum(pred)/sum(y)-1",
       main = paste("Tail calibration, weight", tw))
  abline(h = 0, lty = 2)
  for (i in seq_along(estimators)) {
    g <- z[z$estimator == estimators[i], ]
    lines(g$claim_decile, g$reserve_error_mean, type = "b", pch = 16, col = cols[i])
  }
  legend("bottomleft", legend = estimators, col = cols, pch = 16, lty = 1, cex = 0.75)
}
par(op); invisible(dev.off())

cat(sprintf("\nSaved fremtpl2_tail_sweep_{summary,raw}.csv, fremtpl2_tail_calibration_{summary,raw}.csv, plots | elapsed %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
