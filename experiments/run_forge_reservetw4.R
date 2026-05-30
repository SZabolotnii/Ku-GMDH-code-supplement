#!/usr/bin/env Rscript
# Utah FORGE 58-32 ReserveTW4 cross-geology audit (Paper 1 revision, Phase 1).
#
# Replaces the UCI Gas Turbine CO flagship with the Volve 15/9-F-15 drilling
# benchmark. The reserve window is the OUT-OF-DISTRIBUTION block past the
# gamma-ray formation transition (depth ~2289 m). Pre-transition rows are split
# temporally into train / selection-validation / interval-calibration. The
# ReserveTW4 criterion = reserve-aware external criterion with tail_weight = 4,
# the structure-selection rule that triggered the catastrophic LSE-GMDH failure
# on Gas Turbine; we test whether that failure replicates on drilling data.
#
# Usage:
#   Rscript experiments/run_forge_reservetw4.R [N_SEEDS] [SUBSAMPLE] [B]
#     N_SEEDS   number of repeated-seed runs (default 1 = single-seed audit)
#     SUBSAMPLE pre-transition rows kept per seed (default 0 = all; M2 jitter)
#     B         cumulant-bootstrap resamples for PMM dispatch (default 200)

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS   <- if (length(args) >= 1) as.integer(args[1]) else 1L
SUBSAMPLE <- if (length(args) >= 2) as.integer(args[2]) else 0L
B         <- if (length(args) >= 3) as.integer(args[3]) else 200L
SEED0     <- 75001L

L_MAX <- 3L; F_KEEP <- 6L
TAIL_PROB <- 0.90; TAIL_WEIGHT <- 4; INTERVAL_LEVEL <- 0.95
METHODS <- c("PMM", "LSE", "ridge-LSE", "Huber", "L1")

standardize_by <- function(X, ref) {
  mu <- colMeans(X[ref, , drop = FALSE])
  sig <- apply(X[ref, , drop = FALSE], 2, stats::sd)
  sig[!is.finite(sig) | sig <= 1e-12] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sig, "/")
}

criterion_grid <- function() data.frame(
  label = c("MSE", paste0("ReserveTW", TAIL_WEIGHT)),
  type  = c("MSE", "reserve-aware"),
  reserve_weight = c(0, 1), tail_weight = c(0, TAIL_WEIGHT),
  stringsAsFactors = FALSE)

ctrl_for <- function(method, spec, n_tr, n_val, seed) {
  base <- list(L_max = L_MAX, F = F_KEEP, epsilon = -Inf,
               train_index = seq_len(n_tr), val_index = n_tr + seq_len(n_val),
               criterion = spec$type, reserve_weight = spec$reserve_weight,
               tail_weight = spec$tail_weight, overreserve_weight = 0.25,
               tail_prob = TAIL_PROB, max_iter = 60L)
  switch(method,
    PMM         = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "auto", boot_seed = seed))),
    LSE         = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "LSE"))),
    `ridge-LSE` = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "ridge-LSE", ridge_lambda = 1e-2))),
    Huber       = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "Huber"))),
    L1          = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "L1"))),
    stop("unknown method"))
}

method_share <- function(fit, m) {
  tabs <- lapply(fit$layers, function(z) z$methods)
  num <- sum(vapply(tabs, function(t) { v <- unname(t[m]); if (!length(v) || is.na(v)) 0 else v }, numeric(1)))
  den <- sum(vapply(tabs, sum, numeric(1)))
  if (den <= 0) NA_real_ else num / den
}

max_log10_cond <- function(fit, Xfit) {
  out <- vector("list", length(fit$nodes)); tr <- seq_len(fit$n_train); best <- -Inf
  for (nd in fit$nodes) {
    if (is.null(nd)) next
    if (nd$type == "input") { out[[nd$id]] <- Xfit[, nd$var]; next }
    p1 <- out[[nd$parents[1]]]; p2 <- out[[nd$parents[2]]]
    out[[nd$id]] <- kg2_predict(nd$theta, p1, p2)
    Z <- cbind(1, p1[tr], p2[tr], p1[tr]*p2[tr], p1[tr]^2, p2[tr]^2)
    sv <- svd(Z, nu = 0, nv = 0)$d
    cond <- if (length(sv) && min(sv) > 0) max(sv)/min(sv) else Inf
    best <- max(best, log10(cond))
  }
  best
}

nrmse <- function(y, pred) { s <- stats::sd(y); if (!is.finite(s) || s <= 0) s <- 1; sqrt(mean((pred - y)^2)) / s }

# 95% prediction interval from absolute calibration residuals (symmetric).
interval_cov_width <- function(y_test, pred_test, cal_resid) {
  q <- as.numeric(stats::quantile(abs(cal_resid), INTERVAL_LEVEL, na.rm = TRUE, names = FALSE))
  inside <- y_test >= pred_test - q & y_test <= pred_test + q
  list(coverage = mean(inside), n_test = length(y_test), n_in = sum(inside), width = 2 * q)
}

run_one <- function(D, seed) {
  set.seed(seed)
  T_pre <- D$transition
  pre_idx <- seq_len(T_pre)
  if (SUBSAMPLE > 0 && SUBSAMPLE < T_pre) pre_idx <- sort(sample(pre_idx, SUBSAMPLE))
  T <- length(pre_idx)
  # temporal split of the pre-transition block: 60 / 20 / 20
  n_tr  <- floor(0.60 * T); n_val <- floor(0.20 * T)
  fit_rows_pre <- pre_idx[seq_len(n_tr + n_val)]
  cal_rows     <- pre_idx[(n_tr + n_val + 1L):T]
  test_rows    <- (T_pre + 1L):D$n          # FIXED post-transition reserve window

  Xs <- standardize_by(D$X, fit_rows_pre)
  rows <- list()
  for (i in seq_len(nrow(criterion_grid()))) {
    spec <- criterion_grid()[i, ]
    for (m in METHODS) {
      fit <- gmdh_pmm(Xs[fit_rows_pre, , drop = FALSE], D$y[fit_rows_pre],
                      ctrl_for(m, spec, n_tr, n_val, seed))
      pred_cal  <- predict(fit, Xs[cal_rows, , drop = FALSE])
      pred_test <- predict(fit, Xs[test_rows, , drop = FALSE])
      iv <- interval_cov_width(D$y[test_rows], pred_test, D$y[cal_rows] - pred_cal)
      rows[[length(rows) + 1L]] <- data.frame(
        seed = seed, criterion = spec$label, method = m,
        nrmse_test = nrmse(D$y[test_rows], pred_test),
        nrmse_cal  = nrmse(D$y[cal_rows], pred_cal),
        coverage   = iv$coverage, n_test = iv$n_test, n_in = iv$n_in,
        width = iv$width,
        max_log10_cond = max_log10_cond(fit, Xs[fit_rows_pre, , drop = FALSE]),
        pmm2_share = method_share(fit, "PMM2"), pmm3_share = method_share(fit, "PMM3"),
        best_layer = fit$best$layer, n_train = n_tr, n_val = n_val,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

cat(sprintf("=== Volve ReserveTW4 audit | seeds=%d subsample=%d B=%d ===\n",
            N_SEEDS, SUBSAMPLE, B))
t0 <- Sys.time()
D <- load_forge_drilling(target = "log_rop")
cat(sprintf("FORGE 58-32: n=%d, transition=%d (%.0f%%, %.0f ft), reserve-test=%d rows\n\n",
            D$n, D$transition, 100*D$transition/D$n, D$depth[D$transition], D$n - D$transition))

all <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s) {
  res <- run_one(D, SEED0 + s - 1L)
  cat(sprintf("  seed %d done (%.0fs)\n", SEED0 + s - 1L, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  res
}))

outdir <- "experiments/results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
stem <- if (N_SEEDS > 1) "forge_reservetw4_repeated" else "forge_reservetw4_single"
utils::write.csv(all, file.path(outdir, paste0(stem, "_raw.csv")), row.names = FALSE)

# Aggregate: median + IQR across seeds
agg <- do.call(rbind, by(all, list(all$criterion, all$method), function(g) {
  data.frame(criterion = g$criterion[1], method = g$method[1], n_seed = nrow(g),
    nrmse_med = median(g$nrmse_test), nrmse_q1 = quantile(g$nrmse_test, .25, names = FALSE),
    nrmse_q3 = quantile(g$nrmse_test, .75, names = FALSE), nrmse_max = max(g$nrmse_test),
    cov_med = median(g$coverage), cov_q1 = quantile(g$coverage, .25, names = FALSE),
    cov_q3 = quantile(g$coverage, .75, names = FALSE),
    width_med = median(g$width), cond_med = median(g$max_log10_cond),
    pmm2_med = median(g$pmm2_share, na.rm = TRUE), stringsAsFactors = FALSE) }))
agg <- agg[order(agg$criterion, agg$nrmse_med), ]
utils::write.csv(agg, file.path(outdir, paste0(stem, "_summary.csv")), row.names = FALSE)

cat("\n=== Reserve-test summary (median across seeds) ===\n")
print(agg[, c("criterion","method","n_seed","nrmse_med","nrmse_q1","nrmse_q3","nrmse_max",
              "cov_med","width_med","cond_med","pmm2_med")], row.names = FALSE, digits = 4)
cat(sprintf("\nSaved %s_{raw,summary}.csv | elapsed %.1fs\n", stem,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
