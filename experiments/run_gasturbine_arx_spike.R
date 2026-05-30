#!/usr/bin/env Rscript
# C5 Gas Turbine ARX-like lags and spike-aware validation criterion.
#
# This is a follow-up to run_gasturbine_timeaware.R. It keeps the same
# chronological protocol:
#   internal train: 2011-2013
#   internal validation: 2014
#   external test: 2015
#
# Feature modes:
#   static: current operational inputs only
#   arx:    current operational inputs + target lag1/lag2 + input lag1
#
# Run from the package directory:
#   Rscript experiments/run_gasturbine_arx_spike.R [CO|NOX] [B] [tail_weights]
# tail_weights defaults to 4; pass a comma list for a wider sweep.
# Example:
#   Rscript experiments/run_gasturbine_arx_spike.R CO 30 2,4,8

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
TARGET <- toupper(if (length(args) >= 1) args[1] else "CO")
B <- if (length(args) >= 2) as.integer(args[2]) else 30L
TAIL_WEIGHTS <- if (length(args) >= 3) {
  as.numeric(strsplit(args[3], ",", fixed = TRUE)[[1]])
} else {
  c(4)
}
TAIL_WEIGHTS <- TAIL_WEIGHTS[is.finite(TAIL_WEIGHTS) & TAIL_WEIGHTS >= 0]
if (!(TARGET %in% c("CO", "NOX"))) stop("TARGET must be CO or NOX.")
if (!length(TAIL_WEIGHTS)) stop("tail_weights must contain at least one non-negative number.")

L_MAX <- 3L
F_KEEP <- 6L
PREDICTORS <- c("AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT", "CDP", "TEY")
METHODS <- c("PMM", "LSE", "Huber", "L1")
TAIL_PROB <- 0.90

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

add_arx_lags <- function(d) {
  out <- d
  for (k in 1:2) {
    out[[paste0(TARGET, "_lag", k)]] <- lag_within_year(out[[TARGET]], out$year, k)
  }
  for (nm in PREDICTORS) {
    out[[paste0(nm, "_lag1")]] <- lag_within_year(out[[nm]], out$year, 1L)
  }
  out
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

criterion_grid <- function() {
  rows <- list(data.frame(
    criterion_label = "MSE",
    criterion_type = "MSE",
    reserve_weight = 0,
    tail_weight = 0,
    overreserve_weight = 0.25,
    stringsAsFactors = FALSE
  ))
  for (tw in TAIL_WEIGHTS) {
    rows[[length(rows) + 1L]] <- data.frame(
      criterion_label = paste0("ReserveTW", format(tw, trim = TRUE)),
      criterion_type = "reserve-aware",
      reserve_weight = 1,
      tail_weight = tw,
      overreserve_weight = 0.25,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      criterion_label = paste0("SpikeTW", format(tw, trim = TRUE)),
      criterion_type = "spike-aware",
      reserve_weight = 1,
      tail_weight = tw,
      overreserve_weight = 0.25,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
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
                                                      boot_seed = 73001L))),
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

make_design <- function(d, feature_set) {
  if (feature_set == "static") {
    features <- PREDICTORS
    dd <- d
  } else if (feature_set == "arx") {
    dd <- add_arx_lags(d)
    features <- c(PREDICTORS, paste0(TARGET, "_lag", 1:2),
                  paste0(PREDICTORS, "_lag1"))
  } else {
    stop("Unknown feature_set: ", feature_set)
  }
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
  diag <- bootstrap_cumulant_diag(eps, B = B, robust = TRUE, seed = 73101L)
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
  pred_rows <- list()
  cal_rows <- list()
  for (i in seq_len(nrow(criteria))) {
    spec <- criteria[i, ]
    for (method in METHODS) {
      cat(sprintf("Fitting %-6s | %-10s | %s\n",
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

      pred_rows[[length(pred_rows) + 1L]] <- data.frame(
        target = TARGET,
        feature_set = feature_set,
        method = method,
        criterion_label = spec$criterion_label,
        year = d$year[test_rows],
        row_id = d$row_id[test_rows],
        observed = y_test,
        predicted = pred,
        residual = pred - y_test,
        stringsAsFactors = FALSE
      )
      cal_rows[[length(cal_rows) + 1L]] <- calibration_bins(
        y_test, pred, method, feature_set, spec$criterion_label
      )
    }
  }

  list(
    summary = do.call(rbind, rows),
    predictions = do.call(rbind, pred_rows),
    calibration = do.call(rbind, cal_rows)
  )
}

cat(sprintf("=== Gas Turbine ARX/spike C5 (%s, B=%d, tail_weights=%s) ===\n\n",
            TARGET, B, paste(TAIL_WEIGHTS, collapse = ",")))
t0 <- Sys.time()
root <- find_root()
d <- load_gas_turbine(root)
criteria <- criterion_grid()
static_design <- make_design(d, "static")
arx_design <- make_design(d, "arx")

screen <- rbind(
  residual_screen(static_design, "static"),
  residual_screen(arx_design, "arx")
)
print(screen, row.names = FALSE)

static_criteria <- criteria[criteria$criterion_label == "MSE", , drop = FALSE]
static_run <- run_feature_set(static_design, "static", static_criteria)
arx_run <- run_feature_set(arx_design, "arx", criteria)

summary <- rbind(static_run$summary, arx_run$summary)
summary <- summary[order(summary$feature_set, summary$criterion_label, summary$nrmse), ]
predictions <- rbind(static_run$predictions, arx_run$predictions)
calibration <- rbind(static_run$calibration, arx_run$calibration)

print(summary[, c("feature_set", "criterion_label", "method", "n_features",
                  "nrmse", "mae", "bias", "reserve_error",
                  "top10_mae", "top10_bias", "top10_reserve_error",
                  "top5_mae", "top5_bias", "top5_reserve_error",
                  "layers", "best_layer", "selected_nodes", "pmm2_share")],
      row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
stem <- paste0("gas_turbine_arx_spike_", tolower(TARGET))
utils::write.csv(summary, file.path(outdir, paste0(stem, "_summary.csv")), row.names = FALSE)
utils::write.csv(predictions, file.path(outdir, paste0(stem, "_predictions.csv")), row.names = FALSE)
utils::write.csv(calibration, file.path(outdir, paste0(stem, "_calibration.csv")), row.names = FALSE)
utils::write.csv(screen, file.path(outdir, paste0(stem, "_residual_screen.csv")), row.names = FALSE)

png(file.path(outdir, paste0(stem, ".png")), width = 1280, height = 860, res = 120)
op <- par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))
arx_summary <- summary[summary$feature_set == "arx", ]
labels <- paste(arx_summary$criterion_label, arx_summary$method, sep = "\n")
cols <- ifelse(arx_summary$method == "PMM", "firebrick",
               ifelse(arx_summary$method == "LSE", "grey45",
                      ifelse(arx_summary$method == "Huber", "darkgreen", "purple4")))
barplot(arx_summary$nrmse, names.arg = labels, las = 2, col = cols,
        main = paste(TARGET, "ARX 2015 NRMSE"), ylab = "NRMSE", cex.names = 0.7)
barplot(arx_summary$top10_reserve_error, names.arg = labels, las = 2, col = cols,
        main = "ARX top-10% reserve error", ylab = "sum(pred)/sum(y)-1",
        cex.names = 0.7)
abline(h = 0, lty = 2)
barplot(arx_summary$top5_reserve_error, names.arg = labels, las = 2, col = cols,
        main = "ARX top-5% reserve error", ylab = "sum(pred)/sum(y)-1",
        cex.names = 0.7)
abline(h = 0, lty = 2)
best <- summary[order(summary$nrmse), ][seq_len(min(12L, nrow(summary))), ]
best_labels <- paste(best$feature_set, best$criterion_label, best$method, sep = "\n")
barplot(best$nrmse, names.arg = best_labels, las = 2,
        col = ifelse(best$method == "PMM", "firebrick", "grey60"),
        main = "Best NRMSE configs", ylab = "NRMSE", cex.names = 0.7)
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s_{summary,predictions,calibration,residual_screen}.csv and %s.png | elapsed %.1fs\n",
            stem, stem, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
