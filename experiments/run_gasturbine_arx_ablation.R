#!/usr/bin/env Rscript
# C5 Gas Turbine compact ARX feature ablation.
#
# This follow-up asks whether the Gas Turbine gain comes mostly from target
# memory, input memory, or their combination. It keeps the chronological split:
#   internal train: 2011-2013
#   internal validation: 2014
#   external test: 2015
#
# Feature sets:
#   static:      current operational inputs only
#   target_lags: static + target lag1/lag2
#   input_lags:  static + input lag1
#   full_arx:    static + target lag1/lag2 + input lag1
#
# Run from the package directory:
#   Rscript experiments/run_gasturbine_arx_ablation.R [CO|NOX] [B]

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
FEATURE_SETS <- c("static", "target_lags", "input_lags", "full_arx")
METHODS <- c("PMM", "LSE", "Huber", "L1")
TAIL_PROB <- 0.90
TAIL_WEIGHT <- 4

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

add_lags <- function(d) {
  out <- d
  for (k in 1:2) {
    out[[paste0(TARGET, "_lag", k)]] <- lag_within_year(out[[TARGET]], out$year, k)
  }
  for (nm in PREDICTORS) {
    out[[paste0(nm, "_lag1")]] <- lag_within_year(out[[nm]], out$year, 1L)
  }
  out
}

features_for <- function(feature_set) {
  target_lags <- paste0(TARGET, "_lag", 1:2)
  input_lags <- paste0(PREDICTORS, "_lag1")
  switch(feature_set,
         static = PREDICTORS,
         target_lags = c(PREDICTORS, target_lags),
         input_lags = c(PREDICTORS, input_lags),
         full_arx = c(PREDICTORS, target_lags, input_lags),
         stop("Unknown feature_set: ", feature_set))
}

standardize_by_train <- function(X, train) {
  mu <- colMeans(X[train, , drop = FALSE])
  sig <- apply(X[train, , drop = FALSE], 2, stats::sd)
  sig[!is.finite(sig) | sig <= 1e-12] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sig, "/")
}

make_design <- function(d, feature_set) {
  dd <- add_lags(d)
  features <- features_for(feature_set)
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
                                                      boot_seed = 74001L))),
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

calibration_bins <- function(y, pred, method, feature_set, criterion_label) {
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
      feature_set = feature_set,
      criterion_label = criterion_label,
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

residual_screen <- function(design, feature_set) {
  d <- design$data
  train_rows <- which(d$year %in% 2011:2013)
  val_rows <- which(d$year == 2014)
  X <- design$X
  y <- design$y
  mm_train <- cbind(Intercept = 1, X[train_rows, , drop = FALSE])
  mm_val <- cbind(Intercept = 1, X[val_rows, , drop = FALSE])
  fit <- stats::lm.fit(mm_train, y[train_rows])
  pred_train <- as.numeric(mm_train %*% fit$coefficients)
  pred_val <- as.numeric(mm_val %*% fit$coefficients)
  eps <- y[train_rows] - pred_train
  diag <- bootstrap_cumulant_diag(eps, B = B, robust = TRUE, seed = 74101L)
  dispatch <- dispatch_method(diag)
  cbind(
    data.frame(
      target = TARGET,
      feature_set = feature_set,
      n_features = ncol(X),
      n_train = length(train_rows),
      n_val = length(val_rows),
      dispatch = dispatch,
      gamma3 = diag$gamma3,
      gamma4 = diag$gamma4,
      g2 = diag$g2,
      stringsAsFactors = FALSE
    ),
    metrics(y[val_rows], pred_val, prefix = "val_")
  )
}

run_feature_set <- function(design, feature_set, criteria) {
  d <- design$data
  X <- design$X
  y <- design$y
  train_rows <- which(d$year %in% 2011:2013)
  val_rows <- which(d$year == 2014)
  test_rows <- which(d$year == 2015)
  if (!length(train_rows) || !length(val_rows) || !length(test_rows)) {
    stop("Expected non-empty 2011-2013 train, 2014 validation, and 2015 test partitions.")
  }

  fit_rows <- c(train_rows, val_rows)
  n_train <- length(train_rows)
  n_val <- length(val_rows)
  X_fit <- X[fit_rows, , drop = FALSE]
  y_fit <- y[fit_rows]
  X_test <- X[test_rows, , drop = FALSE]
  y_test <- y[test_rows]

  rows <- list()
  cal_rows <- list()
  for (i in seq_len(nrow(criteria))) {
    spec <- criteria[i, ]
    for (method in METHODS) {
      cat(sprintf("Fitting %-11s | %-10s | %s\n",
                  feature_set, spec$criterion_label, method))
      flush.console()
      fit <- gmdh_pmm(X_fit, y_fit, ctrl_for(method, spec, n_train, n_val))
      pred <- predict(fit, X_test)
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(
          target = TARGET,
          feature_set = feature_set,
          n_features = ncol(X),
          method = method,
          criterion_label = spec$criterion_label,
          criterion_type = spec$criterion_type,
          tail_weight = spec$tail_weight,
          train_years = "2011-2013",
          validation_year = 2014L,
          test_year = 2015L,
          n_train = n_train,
          n_val = n_val,
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
      cal_rows[[length(cal_rows) + 1L]] <- calibration_bins(
        y_test, pred, method, feature_set, spec$criterion_label
      )
    }
  }
  list(summary = do.call(rbind, rows), calibration = do.call(rbind, cal_rows))
}

cat(sprintf("=== Gas Turbine compact ARX ablation (%s, B=%d, ReserveTW%d) ===\n\n",
            TARGET, B, TAIL_WEIGHT))
t0 <- Sys.time()
root <- find_root()
d <- load_gas_turbine(root)
criteria <- criterion_grid()
designs <- stats::setNames(lapply(FEATURE_SETS, function(fs) make_design(d, fs)), FEATURE_SETS)

screen <- do.call(rbind, Map(residual_screen, designs, names(designs)))
print(screen, row.names = FALSE)

runs <- Map(run_feature_set, designs, names(designs), MoreArgs = list(criteria = criteria))
summary <- do.call(rbind, lapply(runs, `[[`, "summary"))
summary <- summary[order(summary$feature_set, summary$criterion_label, summary$nrmse), ]
calibration <- do.call(rbind, lapply(runs, `[[`, "calibration"))

print(summary[, c("feature_set", "criterion_label", "method", "n_features",
                  "nrmse", "mae", "bias", "reserve_error",
                  "top10_mae", "top10_bias", "top10_reserve_error",
                  "top5_mae", "top5_bias", "top5_reserve_error",
                  "layers", "best_layer", "selected_nodes", "pmm2_share")],
      row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
stem <- paste0("gas_turbine_arx_ablation_", tolower(TARGET))
utils::write.csv(summary, file.path(outdir, paste0(stem, "_summary.csv")), row.names = FALSE)
utils::write.csv(calibration, file.path(outdir, paste0(stem, "_calibration.csv")), row.names = FALSE)
utils::write.csv(screen, file.path(outdir, paste0(stem, "_residual_screen.csv")), row.names = FALSE)

png(file.path(outdir, paste0(stem, ".png")), width = 1280, height = 860, res = 120)
op <- par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))
plot_metric <- function(metric, title, ylab) {
  z <- summary[summary$criterion_label == paste0("ReserveTW", TAIL_WEIGHT), ]
  labels <- paste(z$feature_set, z$method, sep = "\n")
  cols <- ifelse(z$method == "PMM", "firebrick",
                 ifelse(z$method == "LSE", "grey45",
                        ifelse(z$method == "Huber", "darkgreen", "purple4")))
  barplot(z[[metric]], names.arg = labels, las = 2, col = cols,
          main = title, ylab = ylab, cex.names = 0.7)
  if (grepl("reserve", metric)) abline(h = 0, lty = 2)
}
plot_metric("nrmse", "ReserveTW4 NRMSE", "NRMSE")
plot_metric("top10_reserve_error", "ReserveTW4 top-10 reserve", "sum(pred)/sum(y)-1")
plot_metric("top5_reserve_error", "ReserveTW4 top-5 reserve", "sum(pred)/sum(y)-1")
pmms <- summary[summary$method == "PMM" & summary$criterion_label == paste0("ReserveTW", TAIL_WEIGHT), ]
barplot(pmms$nrmse, names.arg = pmms$feature_set, las = 2, col = "firebrick",
        main = "PMM feature-set ablation", ylab = "NRMSE", cex.names = 0.8)
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s_{summary,calibration,residual_screen}.csv and %s.png | elapsed %.1fs\n",
            stem, stem, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
