#!/usr/bin/env Rscript
# C5 Gas Turbine time-aware validation and spike calibration.
#
# The random-split smoke test identified CO emissions as a survivor. This
# script uses a chronological split:
#   internal train: 2011-2013
#   internal validation: 2014
#   external test: 2015
#
# Run from the package directory:
#   Rscript experiments/run_gasturbine_timeaware.R [CO|NOX] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
TARGET <- toupper(if (length(args) >= 1) args[1] else "CO")
B <- if (length(args) >= 2) as.integer(args[2]) else 60L
if (!(TARGET %in% c("CO", "NOX"))) stop("TARGET must be CO or NOX.")

L_MAX <- 3L
F_KEEP <- 6L
PREDICTORS <- c("AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT", "CDP", "TEY")

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
  d
}

standardize_by_train <- function(X, train) {
  mu <- colMeans(X[train, , drop = FALSE])
  sig <- apply(X[train, , drop = FALSE], 2, stats::sd)
  sig[!is.finite(sig) | sig <= 1e-12] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sig, "/")
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

ctrl_for <- function(method) {
  base <- list(
    L_max = L_MAX, F = F_KEEP, epsilon = -Inf,
    train_index = seq_len(.n_train),
    val_index = .n_train + seq_len(.n_val),
    max_iter = 60L
  )
  switch(method,
         PMM = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "auto",
                                                      boot_seed = 72001L))),
         LSE = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "LSE"))),
         `ridge-LSE` = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "ridge-LSE"))),
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

calibration_bins <- function(y, pred, method) {
  n <- length(y)
  ord <- order(y)
  bin <- integer(n)
  bin[ord] <- pmin(10L, ceiling(seq_len(n) * 10 / n))
  rows <- lapply(seq_len(10L), function(b) {
    idx <- bin == b
    observed_sum <- sum(y[idx])
    predicted_sum <- sum(pred[idx])
    data.frame(
      method = method,
      emission_decile = b,
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

cat(sprintf("=== Gas Turbine time-aware C5 (%s, B=%d) ===\n\n", TARGET, B))
t0 <- Sys.time()
root <- find_root()
d <- load_gas_turbine(root)

train_rows <- which(d$year %in% 2011:2013)
val_rows <- which(d$year == 2014)
test_rows <- which(d$year == 2015)
if (!length(train_rows) || !length(val_rows) || !length(test_rows)) {
  stop("Expected non-empty 2011-2013 train, 2014 validation, and 2015 test partitions.")
}

X <- as.matrix(d[, PREDICTORS])
storage.mode(X) <- "double"
varied <- apply(X[train_rows, , drop = FALSE], 2, stats::sd) > 1e-12
X <- X[, varied, drop = FALSE]
X <- standardize_by_train(X, train_rows)
y <- d[[TARGET]]

fit_rows <- c(train_rows, val_rows)
.n_train <- length(train_rows)
.n_val <- length(val_rows)
X_fit <- X[fit_rows, , drop = FALSE]
y_fit <- y[fit_rows]
X_test <- X[test_rows, , drop = FALSE]
y_test <- y[test_rows]

methods <- c("PMM", "LSE", "ridge-LSE", "Huber", "L1")
rows <- list()
pred_rows <- list()
cal_rows <- list()
for (method in methods) {
  fit <- gmdh_pmm(X_fit, y_fit, ctrl_for(method))
  pred <- predict(fit, X_test)

  rows[[length(rows) + 1L]] <- cbind(
    data.frame(
      target = TARGET,
      method = method,
      train_years = "2011-2013",
      validation_year = 2014L,
      test_year = 2015L,
      n_train = .n_train,
      n_val = .n_val,
      n_test = length(test_rows),
      layers = length(fit$layers),
      best_layer = fit$best$layer,
      selected_nodes = selected_node_count(fit),
      pmm2_share = method_share(fit, "PMM2"),
      pmm3_share = method_share(fit, "PMM3"),
      stringsAsFactors = FALSE
    ),
    metrics(y_test, pred),
    spike_metrics(y_test, pred, 0.90, "top10_"),
    spike_metrics(y_test, pred, 0.95, "top5_")
  )

  pred_rows[[length(pred_rows) + 1L]] <- data.frame(
    target = TARGET,
    method = method,
    year = d$year[test_rows],
    observed = y_test,
    predicted = pred,
    residual = pred - y_test,
    stringsAsFactors = FALSE
  )
  cal_rows[[length(cal_rows) + 1L]] <- calibration_bins(y_test, pred, method)
}

summary <- do.call(rbind, rows)
summary <- summary[order(summary$nrmse), ]
predictions <- do.call(rbind, pred_rows)
calibration <- do.call(rbind, cal_rows)

print(summary[, c("method", "nrmse", "mae", "bias", "reserve_error",
                  "top10_mae", "top10_bias", "top10_reserve_error",
                  "top5_mae", "top5_bias", "top5_reserve_error",
                  "layers", "best_layer", "selected_nodes", "pmm2_share")],
      row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
stem <- paste0("gas_turbine_timeaware_", tolower(TARGET))
utils::write.csv(summary, file.path(outdir, paste0(stem, "_summary.csv")), row.names = FALSE)
utils::write.csv(predictions, file.path(outdir, paste0(stem, "_predictions.csv")), row.names = FALSE)
utils::write.csv(calibration, file.path(outdir, paste0(stem, "_calibration.csv")), row.names = FALSE)

png(file.path(outdir, paste0(stem, ".png")), width = 1120, height = 760, res = 115)
op <- par(mfrow = c(2, 2), mar = c(6, 4, 3, 1))
cols <- ifelse(summary$method == "PMM", "firebrick",
               ifelse(summary$method == "LSE", "grey45",
                      ifelse(summary$method == "ridge-LSE", "steelblue",
                             ifelse(summary$method == "Huber", "darkgreen", "purple4"))))
barplot(summary$nrmse, names.arg = summary$method, las = 2, col = cols,
        main = paste(TARGET, "2015 NRMSE"), ylab = "NRMSE")
barplot(summary$mae, names.arg = summary$method, las = 2, col = cols,
        main = "MAE", ylab = "MAE")
barplot(summary$top10_reserve_error, names.arg = summary$method, las = 2, col = cols,
        main = "Top-10% reserve error", ylab = "sum(pred)/sum(y)-1")
abline(h = 0, lty = 2)
barplot(summary$selected_nodes, names.arg = summary$method, las = 2, col = cols,
        main = "Selected nodes", ylab = "count")
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s_{summary,predictions,calibration}.csv and %s.png | elapsed %.1fs\n",
            stem, stem, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
