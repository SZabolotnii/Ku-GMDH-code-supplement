#!/usr/bin/env Rscript
# C5 Gas Turbine prediction-interval coverage reading.
#
# Real-world CO data do not have known true coefficients, so this is not the
# synthetic H1 coefficient-coverage experiment. Instead, we calibrate prediction
# intervals on 2014 validation residuals and test them on the external 2015
# holdout. The intended reading is practical: do PMM-selected dynamic soft
# sensors keep emission spikes inside calibrated intervals, and how wide are
# those intervals relative to LSE-family baselines?
#
# Run from the package directory:
#   Rscript experiments/run_gasturbine_interval_coverage.R [CO|NOX] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
TARGET <- toupper(if (length(args) >= 1) args[1] else "CO")
B <- if (length(args) >= 2) as.integer(args[2]) else 30L
if (!(TARGET %in% c("CO", "NOX"))) stop("TARGET must be CO or NOX.")

L_MAX <- 3L
F_KEEP <- 6L
PREDICTORS <- c("AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT", "CDP", "TEY")
METHODS <- c("PMM", "LSE", "Huber", "L1")
TAIL_PROB <- 0.90
TAIL_WEIGHT <- 4
FEATURE_SET <- "target_lags"
INTERVAL_LEVEL <- 0.95

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

load_gas_turbine <- function(root) {
  path <- file.path(root, "shared", "datasets", "external", "processed",
                    "gas_turbine_2011_2015.csv")
  if (!file.exists(path)) {
    stop("Missing prepared Gas Turbine data: ", path,
         "\nRun: Rscript shared/datasets/prepare_gasturbine.R")
  }
  d <- stats::na.omit(utils::read.csv(path, stringsAsFactors = FALSE))
  required <- c("year", TARGET, PREDICTORS)
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  d$row_id <- seq_len(nrow(d))
  d
}

lag_within_year <- function(x, year, k) {
  out <- rep(NA_real_, length(x))
  for (yy in sort(unique(year))) {
    idx <- which(year == yy)
    if (length(idx) > k) out[idx[(k + 1L):length(idx)]] <- x[idx[seq_len(length(idx) - k)]]
  }
  out
}

add_target_lags <- function(d) {
  out <- d
  for (k in 1:2) {
    out[[paste0(TARGET, "_lag", k)]] <- lag_within_year(out[[TARGET]], out$year, k)
  }
  out
}

standardize_by_train <- function(X, train) {
  mu <- colMeans(X[train, , drop = FALSE])
  sig <- apply(X[train, , drop = FALSE], 2, stats::sd)
  sig[!is.finite(sig) | sig <= 1e-12] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sig, "/")
}

make_design <- function(d) {
  dd <- add_target_lags(d)
  features <- c(PREDICTORS, paste0(TARGET, "_lag", 1:2))
  keep <- stats::complete.cases(dd[, c("year", TARGET, features)])
  dd <- dd[keep, , drop = FALSE]
  X <- as.matrix(dd[, features])
  storage.mode(X) <- "double"
  train_rows <- which(dd$year %in% 2011:2013)
  varied <- apply(X[train_rows, , drop = FALSE], 2, stats::sd) > 1e-12
  X <- X[, varied, drop = FALSE]
  X <- standardize_by_train(X, train_rows)
  list(data = dd, X = X, y = dd[[TARGET]], features = colnames(X))
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

selected_node_count <- function(fit) {
  sum(vapply(fit$layers, function(z) sum(z$methods), numeric(1)))
}

criterion_grid <- function() {
  data.frame(
    criterion_label = c("MSE", paste0("ReserveTW", TAIL_WEIGHT)),
    criterion_type = c("MSE", "reserve-aware"),
    reserve_weight = c(0, 1),
    tail_weight = c(0, TAIL_WEIGHT),
    overreserve_weight = c(0.25, 0.25),
    stringsAsFactors = FALSE
  )
}

ctrl_for <- function(method, spec, n_train, n_val) {
  base <- list(
    L_max = L_MAX, F = F_KEEP, epsilon = -Inf,
    train_index = seq_len(n_train),
    val_index = n_train + seq_len(n_val),
    criterion = spec$criterion_type,
    reserve_weight = spec$reserve_weight,
    tail_weight = spec$tail_weight,
    overreserve_weight = spec$overreserve_weight,
    tail_prob = TAIL_PROB,
    max_iter = 60L
  )
  switch(method,
         PMM = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "auto",
                                                      boot_seed = 75001L))),
         LSE = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "LSE"))),
         Huber = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "Huber"))),
         L1 = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "L1"))),
         stop("Unknown method: ", method))
}

metrics <- function(y, pred, prefix = "") {
  e <- pred - y
  y_scale <- stats::sd(y)
  if (!is.finite(y_scale) || y_scale <= 1e-12) y_scale <- 1
  out <- data.frame(
    nrmse = sqrt(mean(e^2)) / y_scale,
    mae = mean(abs(e)),
    bias = mean(e),
    under_rate = mean(pred < y),
    reserve_error = sum(pred) / sum(y) - 1,
    stringsAsFactors = FALSE
  )
  names(out) <- paste0(prefix, names(out))
  out
}

spike_metrics <- function(y, pred, prob, prefix) {
  threshold <- as.numeric(stats::quantile(y, prob, names = FALSE, na.rm = TRUE))
  idx <- y >= threshold
  head <- data.frame(threshold = threshold, n = sum(idx), stringsAsFactors = FALSE)
  names(head) <- paste0(prefix, names(head))
  cbind(head, metrics(y[idx], pred[idx], prefix = prefix))
}

interval_specs <- function() {
  alpha <- 1 - INTERVAL_LEVEL
  data.frame(
    interval_type = c("equal_tail_95", "abs_residual_95", "upper_95"),
    lower_prob = c(alpha / 2, NA_real_, NA_real_),
    upper_prob = c(1 - alpha / 2, NA_real_, INTERVAL_LEVEL),
    abs_prob = c(NA_real_, INTERVAL_LEVEL, NA_real_),
    one_sided = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

interval_offsets <- function(cal_resid, spec) {
  if (identical(spec$interval_type, "abs_residual_95")) {
    q <- as.numeric(stats::quantile(abs(cal_resid), spec$abs_prob,
                                    na.rm = TRUE, names = FALSE))
    return(list(lower_offset = -q, upper_offset = q, width = 2 * q))
  }
  if (isTRUE(spec$one_sided)) {
    q <- as.numeric(stats::quantile(cal_resid, spec$upper_prob,
                                    na.rm = TRUE, names = FALSE))
    return(list(lower_offset = -Inf, upper_offset = q, width = NA_real_))
  }
  lo <- as.numeric(stats::quantile(cal_resid, spec$lower_prob,
                                   na.rm = TRUE, names = FALSE))
  hi <- as.numeric(stats::quantile(cal_resid, spec$upper_prob,
                                   na.rm = TRUE, names = FALSE))
  list(lower_offset = lo, upper_offset = hi, width = hi - lo)
}

score_interval <- function(y, pred, lower, upper, prefix = "") {
  in_interval <- y >= lower & y <= upper
  lower_miss <- y < lower
  upper_miss <- y > upper
  width <- upper - lower
  width[!is.finite(width)] <- NA_real_
  out <- data.frame(
    coverage = mean(in_interval),
    lower_miss_rate = mean(lower_miss),
    upper_miss_rate = mean(upper_miss),
    mean_width = mean(width, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  names(out) <- paste0(prefix, names(out))
  out
}

tail_interval_scores <- function(y, pred, lower, upper) {
  top10 <- y >= as.numeric(stats::quantile(y, 0.90, na.rm = TRUE, names = FALSE))
  top5 <- y >= as.numeric(stats::quantile(y, 0.95, na.rm = TRUE, names = FALSE))
  cbind(score_interval(y[top10], pred[top10], lower[top10], upper[top10], "top10_"),
        score_interval(y[top5], pred[top5], lower[top5], upper[top5], "top5_"))
}

coverage_by_decile <- function(y, pred, lower, upper, method, criterion_label, interval_type) {
  n <- length(y)
  ord <- order(y)
  bin <- integer(n)
  bin[ord] <- pmin(10L, ceiling(seq_len(n) * 10 / n))
  rows <- lapply(seq_len(10L), function(b) {
    idx <- bin == b
    width <- upper[idx] - lower[idx]
    width[!is.finite(width)] <- NA_real_
    data.frame(
      method = method,
      criterion_label = criterion_label,
      interval_type = interval_type,
      emission_decile = b,
      n = sum(idx),
      observed_mean = mean(y[idx]),
      predicted_mean = mean(pred[idx]),
      coverage = mean(y[idx] >= lower[idx] & y[idx] <= upper[idx]),
      upper_miss_rate = mean(y[idx] > upper[idx]),
      mean_width = mean(width, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fit_and_score <- function(design, method, spec) {
  d <- design$data
  X <- design$X
  y <- design$y
  train_rows <- which(d$year %in% 2011:2013)
  val_rows <- which(d$year == 2014)
  test_rows <- which(d$year == 2015)
  fit_rows <- c(train_rows, val_rows)
  n_train <- length(train_rows)
  n_val <- length(val_rows)

  cat(sprintf("Fitting %-10s | %s\n", spec$criterion_label, method))
  flush.console()
  fit <- gmdh_pmm(X[fit_rows, , drop = FALSE], y[fit_rows],
                  ctrl_for(method, spec, n_train, n_val))

  pred_val <- predict(fit, X[val_rows, , drop = FALSE])
  pred_test <- predict(fit, X[test_rows, , drop = FALSE])
  cal_resid <- y[val_rows] - pred_val

  interval_rows <- list()
  decile_rows <- list()
  for (i in seq_len(nrow(interval_specs()))) {
    ispec <- interval_specs()[i, ]
    off <- interval_offsets(cal_resid, ispec)
    lower_test <- pred_test + off$lower_offset
    upper_test <- pred_test + off$upper_offset
    lower_val <- pred_val + off$lower_offset
    upper_val <- pred_val + off$upper_offset

    interval_rows[[length(interval_rows) + 1L]] <- cbind(
      data.frame(
        target = TARGET,
        feature_set = FEATURE_SET,
        n_features = ncol(X),
        method = method,
        criterion_label = spec$criterion_label,
        criterion_type = spec$criterion_type,
        interval_type = ispec$interval_type,
        level = INTERVAL_LEVEL,
        lower_offset = off$lower_offset,
        upper_offset = off$upper_offset,
        interval_width = off$width,
        train_years = "2011-2013",
        calibration_year = 2014L,
        test_year = 2015L,
        n_train = n_train,
        n_cal = n_val,
        n_test = length(test_rows),
        layers = length(fit$layers),
        best_layer = fit$best$layer,
        selected_nodes = selected_node_count(fit),
        pmm2_share = method_share(fit, "PMM2"),
        pmm3_share = method_share(fit, "PMM3"),
        stringsAsFactors = FALSE
      ),
      metrics(y[test_rows], pred_test),
      spike_metrics(y[test_rows], pred_test, 0.90, "top10_"),
      spike_metrics(y[test_rows], pred_test, 0.95, "top5_"),
      score_interval(y[val_rows], pred_val, lower_val, upper_val, "cal_"),
      score_interval(y[test_rows], pred_test, lower_test, upper_test, "test_"),
      tail_interval_scores(y[test_rows], pred_test, lower_test, upper_test)
    )

    decile_rows[[length(decile_rows) + 1L]] <- coverage_by_decile(
      y[test_rows], pred_test, lower_test, upper_test,
      method, spec$criterion_label, ispec$interval_type
    )
  }

  list(summary = do.call(rbind, interval_rows),
       deciles = do.call(rbind, decile_rows))
}

cat(sprintf("=== Gas Turbine interval coverage (%s, B=%d, feature_set=%s) ===\n\n",
            TARGET, B, FEATURE_SET))
t0 <- Sys.time()
root <- find_root()
d <- load_gas_turbine(root)
design <- make_design(d)
criteria <- criterion_grid()

runs <- list()
for (i in seq_len(nrow(criteria))) {
  spec <- criteria[i, ]
  for (method in METHODS) {
    runs[[length(runs) + 1L]] <- fit_and_score(design, method, spec)
  }
}

summary <- do.call(rbind, lapply(runs, `[[`, "summary"))
summary <- summary[order(summary$criterion_label, summary$interval_type,
                         summary$test_upper_miss_rate, summary$nrmse), ]
deciles <- do.call(rbind, lapply(runs, `[[`, "deciles"))

print(summary[, c("criterion_label", "method", "interval_type",
                  "nrmse", "mae", "reserve_error",
                  "top10_reserve_error", "top5_reserve_error",
                  "cal_coverage", "test_coverage", "test_upper_miss_rate",
                  "test_mean_width", "top10_coverage", "top10_upper_miss_rate",
                  "top5_coverage", "top5_upper_miss_rate",
                  "pmm2_share")],
      row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
stem <- paste0("gas_turbine_interval_coverage_", tolower(TARGET))
utils::write.csv(summary, file.path(outdir, paste0(stem, "_summary.csv")), row.names = FALSE)
utils::write.csv(deciles, file.path(outdir, paste0(stem, "_deciles.csv")), row.names = FALSE)

png(file.path(outdir, paste0(stem, ".png")), width = 1280, height = 860, res = 120)
op <- par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))
plot_block <- function(interval_type, metric, title, ylab) {
  z <- summary[summary$criterion_label == paste0("ReserveTW", TAIL_WEIGHT) &
                 summary$interval_type == interval_type, ]
  labels <- z$method
  cols <- ifelse(z$method == "PMM", "firebrick",
                 ifelse(z$method == "LSE", "grey45",
                        ifelse(z$method == "Huber", "darkgreen", "purple4")))
  barplot(z[[metric]], names.arg = labels, las = 2, col = cols,
          main = title, ylab = ylab, cex.names = 0.8)
  if (grepl("coverage", metric)) abline(h = INTERVAL_LEVEL, lty = 2, col = "grey30")
}
plot_block("equal_tail_95", "test_coverage", "Equal-tail 95% test coverage", "coverage")
plot_block("equal_tail_95", "test_mean_width", "Equal-tail 95% mean width", "width")
plot_block("upper_95", "test_coverage", "Upper 95% test coverage", "coverage")
plot_block("upper_95", "top10_upper_miss_rate", "Upper 95% top-10 miss rate", "miss rate")
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s_{summary,deciles}.csv and %s.png | elapsed %.1fs\n",
            stem, stem, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
