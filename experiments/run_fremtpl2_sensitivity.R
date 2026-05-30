#!/usr/bin/env Rscript
# C5 freMTPL2 sensitivity experiment.
#
# This script follows up the external freMTPL2 candidate-stage result by
# separating three questions:
#   1. default PMM auto-dispatch;
#   2. relaxed PMM2 dispatch under extreme raw severity asymmetry;
#   3. forced PMM2 in every KG-2 partial model.
#
# It reports prediction metrics and actuarial severity diagnostics: aggregate
# reserve error and tail reserve error. The goal is not to make PMM win every
# generic metric, but to test whether PMM changes the under-reserving behavior
# seen with robust losses.
#
# Run from the package directory:
#   Rscript experiments/run_fremtpl2_sensitivity.R [R] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
R_REP <- if (length(args) >= 1) as.integer(args[1]) else 10L
B <- if (length(args) >= 2) as.integer(args[2]) else 40L
TEST_FRAC <- 0.25
L_MAX <- 3L
F_KEEP <- 6L

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
    label = "freMTPL2 raw claim severity",
    y = stats::model.response(mf),
    X = model_x(form, mf)
  )
}

ctrl_for <- function(arm, seed) {
  base <- list(L_max = L_MAX, F = F_KEEP, epsilon = -Inf, seed = seed,
               max_iter = 60L)
  switch(arm,
         `PMM-auto` = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "auto"))),
         `PMM-relaxed` = do.call(gmdh_pmm_control, c(base, list(
           B = B, force_method = "auto", alpha = 0.10,
           skew_min = 0.10, skew_strong = 0.50, g2_threshold = 0.99
         ))),
         `PMM2-forced` = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "PMM2"))),
         LSE = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "LSE"))),
         `ridge-LSE` = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "ridge-LSE"))),
         Huber = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "Huber"))),
         L1 = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "L1"))),
         stop("Unknown arm: ", arm))
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
  tail <- y >= as.numeric(stats::quantile(y, 0.90, na.rm = TRUE))
  data.frame(
    nrmse = sqrt(mean(e^2)) / y_scale,
    mae = mean(abs(e)),
    bias = mean(e),
    reserve_error = sum(pred) / sum(y) - 1,
    under_rate = mean(pred < y),
    tail_mae = mean(abs(e[tail])),
    tail_bias = mean(e[tail]),
    tail_reserve_error = sum(pred[tail]) / sum(y[tail]) - 1,
    tail_under_rate = mean(pred[tail] < y[tail]),
    stringsAsFactors = FALSE
  )
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
  arms <- c("PMM-auto", "PMM-relaxed", "PMM2-forced", "LSE", "ridge-LSE", "Huber", "L1")

  rows <- list()
  for (arm in arms) {
    fit <- try(gmdh_pmm(xs$X_train, y_train, ctrl_for(arm, seed)), silent = TRUE)
    if (inherits(fit, "try-error")) {
      rows[[length(rows) + 1L]] <- data.frame(
        candidate = dat$id, rep = rep_id, arm = arm, ok = FALSE,
        nrmse = NA_real_, mae = NA_real_, bias = NA_real_,
        reserve_error = NA_real_, under_rate = NA_real_,
        tail_mae = NA_real_, tail_bias = NA_real_,
        tail_reserve_error = NA_real_, tail_under_rate = NA_real_,
        pmm2_share = NA_real_, pmm3_share = NA_real_, layers = NA_integer_,
        stringsAsFactors = FALSE)
      next
    }
    pred <- predict(fit, xs$X_test)
    m <- metrics(y_test, pred)
    rows[[length(rows) + 1L]] <- data.frame(
      candidate = dat$id, rep = rep_id, arm = arm, ok = TRUE,
      m,
      pmm2_share = method_share(fit, "PMM2"),
      pmm3_share = method_share(fit, "PMM3"),
      layers = length(fit$layers),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$nrmse_winner <- out$arm == out$arm[which.min(out$nrmse)]
  out$reserve_winner <- abs(out$reserve_error) == min(abs(out$reserve_error), na.rm = TRUE)
  out
}

summarize <- function(raw) {
  parts <- split(raw[raw$ok, ], raw$arm[raw$ok], drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(g) {
    data.frame(
      arm = g$arm[1],
      nrmse_mean = round(mean(g$nrmse), 3),
      nrmse_sd = round(stats::sd(g$nrmse), 3),
      mae_mean = round(mean(g$mae), 1),
      bias_mean = round(mean(g$bias), 1),
      reserve_error_mean = round(mean(g$reserve_error), 3),
      under_rate_mean = round(mean(g$under_rate), 3),
      tail_bias_mean = round(mean(g$tail_bias), 1),
      tail_reserve_error_mean = round(mean(g$tail_reserve_error), 3),
      tail_under_rate_mean = round(mean(g$tail_under_rate), 3),
      nrmse_win_rate = round(mean(g$nrmse_winner), 3),
      reserve_win_rate = round(mean(g$reserve_winner), 3),
      pmm2_share = round(mean(g$pmm2_share), 3),
      layers_mean = round(mean(g$layers), 2),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$nrmse_mean), ]
}

cat(sprintf("=== freMTPL2 sensitivity (R=%d, B=%d) ===\n\n", R_REP, B))
t0 <- Sys.time()
root <- find_root()
dat <- load_fremtpl2_raw(root)
raw <- do.call(rbind, lapply(seq_len(R_REP), function(r) {
  one_rep(dat, r, 92000L + r)
}))
res <- summarize(raw)
print(res, row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(raw, file.path(outdir, "fremtpl2_sensitivity_raw.csv"), row.names = FALSE)
utils::write.csv(res, file.path(outdir, "fremtpl2_sensitivity_summary.csv"), row.names = FALSE)

png(file.path(outdir, "fremtpl2_sensitivity.png"), width = 980, height = 720, res = 115)
op <- par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))
cols <- ifelse(grepl("PMM", res$arm), "firebrick",
               ifelse(res$arm %in% c("Huber", "L1"), "darkgreen", "grey50"))
barplot(res$nrmse_mean, names.arg = res$arm, las = 2, col = cols,
        main = "Test NRMSE", ylab = "mean")
barplot(res$mae_mean, names.arg = res$arm, las = 2, col = cols,
        main = "MAE", ylab = "mean")
barplot(res$reserve_error_mean, names.arg = res$arm, las = 2, col = cols,
        main = "Aggregate reserve error", ylab = "sum(pred)/sum(y)-1")
abline(h = 0, lty = 2)
barplot(res$tail_reserve_error_mean, names.arg = res$arm, las = 2, col = cols,
        main = "Top-decile reserve error", ylab = "tail sum(pred)/sum(y)-1")
abline(h = 0, lty = 2)
par(op); invisible(dev.off())

cat(sprintf("\nSaved fremtpl2_sensitivity_summary.csv, fremtpl2_sensitivity_raw.csv, fremtpl2_sensitivity.png | elapsed %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
