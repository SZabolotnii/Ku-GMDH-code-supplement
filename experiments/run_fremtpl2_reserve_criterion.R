#!/usr/bin/env Rscript
# C5 freMTPL2 reserve-aware external criterion experiment.
#
# This experiment changes the GMDH external selection criterion, not the final
# reporting metrics. It compares classical MSE selection against two
# reserve-aware variants that penalize aggregate and top-decile under-reserve.
#
# Run from the package directory:
#   Rscript experiments/run_fremtpl2_reserve_criterion.R [R] [B]

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

criterion_specs <- list(
  MSE = list(criterion = "MSE", reserve_weight = 0, tail_weight = 0,
             overreserve_weight = 0.25),
  Reserve = list(criterion = "reserve-aware", reserve_weight = 1, tail_weight = 1,
                 overreserve_weight = 0.25),
  TailReserve = list(criterion = "reserve-aware", reserve_weight = 1, tail_weight = 4,
                     overreserve_weight = 0.25)
)

ctrl_for <- function(estimator, criterion_name, seed) {
  spec <- criterion_specs[[criterion_name]]
  base <- list(
    L_max = L_MAX, F = F_KEEP, epsilon = -Inf, seed = seed,
    criterion = spec$criterion,
    reserve_weight = spec$reserve_weight,
    tail_weight = spec$tail_weight,
    overreserve_weight = spec$overreserve_weight,
    tail_prob = 0.9,
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
  tail <- y >= as.numeric(stats::quantile(y, 0.90, na.rm = TRUE))
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
  criteria <- names(criterion_specs)

  rows <- list()
  for (criterion_name in criteria) {
    for (estimator in estimators) {
      fit <- try(gmdh_pmm(xs$X_train, y_train,
                          ctrl_for(estimator, criterion_name, seed)),
                 silent = TRUE)
      if (inherits(fit, "try-error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          rep = rep_id, criterion = criterion_name, estimator = estimator, ok = FALSE,
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
        rep = rep_id, criterion = criterion_name, estimator = estimator, ok = TRUE,
        m,
        pmm2_share = method_share(fit, "PMM2"),
        layers = length(fit$layers),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$nrmse_winner <- ave(out$nrmse, out$rep, out$criterion,
                          FUN = function(z) z == min(z, na.rm = TRUE)) > 0
  out$reserve_winner <- ave(out$abs_reserve_error, out$rep, out$criterion,
                            FUN = function(z) z == min(z, na.rm = TRUE)) > 0
  out$tail_reserve_winner <- ave(out$abs_tail_reserve_error, out$rep, out$criterion,
                                 FUN = function(z) z == min(z, na.rm = TRUE)) > 0
  out
}

summarize <- function(raw) {
  parts <- split(raw[raw$ok, ], list(raw$criterion[raw$ok], raw$estimator[raw$ok]), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(g) {
    data.frame(
      criterion = g$criterion[1],
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
  out[order(out$criterion, out$nrmse_mean), ]
}

cat(sprintf("=== freMTPL2 reserve-aware criterion (R=%d, B=%d) ===\n\n", R_REP, B))
t0 <- Sys.time()
root <- find_root()
dat <- load_fremtpl2_raw(root)
raw <- do.call(rbind, lapply(seq_len(R_REP), function(r) {
  one_rep(dat, r, 104000L + r)
}))
res <- summarize(raw)
print(res, row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(raw, file.path(outdir, "fremtpl2_reserve_criterion_raw.csv"), row.names = FALSE)
utils::write.csv(res, file.path(outdir, "fremtpl2_reserve_criterion_summary.csv"), row.names = FALSE)

png(file.path(outdir, "fremtpl2_reserve_criterion.png"), width = 1080, height = 720, res = 115)
op <- par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))
estimator_order <- c("PMM-auto", "PMM-relaxed", "LSE", "Huber", "L1")
criterion_order <- names(criterion_specs)
metric_matrix <- function(metric) {
  sapply(criterion_order, function(cr) {
    z <- res[res$criterion == cr, ]
    vals <- stats::setNames(z[[metric]], z$estimator)
    vals[estimator_order]
  })
}
for (metric in c("nrmse_mean", "abs_reserve_error_mean", "abs_tail_reserve_error_mean", "mae_mean")) {
  mat <- metric_matrix(metric)
  rownames(mat) <- estimator_order
  barplot(t(mat), beside = TRUE, las = 2, ylab = metric,
          main = metric, col = c("grey55", "firebrick", "steelblue"))
  legend("topright", legend = colnames(mat),
         fill = c("grey55", "firebrick", "steelblue"), cex = 0.75)
}
par(op); invisible(dev.off())

cat(sprintf("\nSaved fremtpl2_reserve_criterion_summary.csv, fremtpl2_reserve_criterion_raw.csv, fremtpl2_reserve_criterion.png | elapsed %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
