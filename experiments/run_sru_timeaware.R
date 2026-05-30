#!/usr/bin/env Rscript
# C5 SRU soft-sensor chronological validation.
#
# SRU is a single normalized MIMO process sequence. This script avoids random
# leakage by using the first 60% of rows for internal training, the next 20%
# for GMDH validation, and the final 20% as external holdout.
#
# Run from the package directory:
#   Rscript experiments/run_sru_timeaware.R [y1|y2] [static|dynamic|both] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
TARGET <- tolower(if (length(args) >= 1) args[1] else "y1")
FEATURE_MODE <- tolower(if (length(args) >= 2) args[2] else "both")
B <- if (length(args) >= 3) as.integer(args[3]) else 30L
if (!(TARGET %in% c("y1", "y2"))) stop("TARGET must be y1 or y2.")
if (!(FEATURE_MODE %in% c("static", "dynamic", "both"))) {
  stop("FEATURE_MODE must be static, dynamic, or both.")
}

L_MAX <- 3L
F_KEEP <- 6L
INPUTS <- paste0("u", 1:5)
METHODS <- c("PMM", "LSE", "ridge-LSE", "Huber", "L1")

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

load_sru <- function(root, feature_set) {
  file_name <- if (feature_set == "static") "sru.csv" else "sru_lagged.csv"
  path <- file.path(root, "shared", "datasets", "external", "processed", file_name)
  if (!file.exists(path)) {
    stop("Missing prepared SRU data: ", path,
         "\nRun: Rscript shared/datasets/prepare_sru.R")
  }
  d <- stats::na.omit(utils::read.csv(path, stringsAsFactors = FALSE))
  required <- c("time", INPUTS, "y1", "y2")
  if (feature_set == "dynamic") {
    required <- c(required, "y1_lag1", "y1_lag2", "y2_lag1", "y2_lag2",
                  paste0(INPUTS, "_lag1"))
  }
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  d[order(d$time), ]
}

features_for <- function(feature_set) {
  if (feature_set == "static") return(INPUTS)
  c(INPUTS, "y1_lag1", "y1_lag2", "y2_lag1", "y2_lag2", paste0(INPUTS, "_lag1"))
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

ctrl_for <- function(method, n_train, n_val) {
  base <- list(
    L_max = L_MAX, F = F_KEEP, epsilon = -Inf,
    train_index = seq_len(n_train),
    val_index = n_train + seq_len(n_val),
    max_iter = 60L
  )
  switch(method,
         PMM = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "auto",
                                                      boot_seed = 76001L))),
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

tail_metrics <- function(y, pred, prob, prefix) {
  threshold <- as.numeric(stats::quantile(y, prob, names = FALSE, na.rm = TRUE))
  idx <- y >= threshold
  head <- data.frame(threshold = threshold, n = sum(idx), stringsAsFactors = FALSE)
  names(head) <- paste0(prefix, names(head))
  cbind(head, metrics(y[idx], pred[idx], prefix = prefix))
}

residual_screen <- function(X, y, train_rows, val_rows, feature_set) {
  mm_train <- cbind(Intercept = 1, X[train_rows, , drop = FALSE])
  mm_val <- cbind(Intercept = 1, X[val_rows, , drop = FALSE])
  fit <- stats::lm.fit(mm_train, y[train_rows])
  pred_train <- as.numeric(mm_train %*% fit$coefficients)
  pred_val <- as.numeric(mm_val %*% fit$coefficients)
  eps <- y[train_rows] - pred_train
  diag <- bootstrap_cumulant_diag(eps, B = B, robust = TRUE, seed = 76101L)
  cbind(
    data.frame(
      target = TARGET,
      feature_set = feature_set,
      n_features = ncol(X),
      n_train = length(train_rows),
      n_val = length(val_rows),
      dispatch = dispatch_method(diag),
      gamma3 = diag$gamma3,
      gamma4 = diag$gamma4,
      g2 = diag$g2,
      g3 = diag$g3,
      stringsAsFactors = FALSE
    ),
    metrics(y[val_rows], pred_val, prefix = "val_")
  )
}

run_feature_set <- function(root, feature_set) {
  d <- load_sru(root, feature_set)
  features <- features_for(feature_set)
  X <- as.matrix(d[, features])
  storage.mode(X) <- "double"
  y <- d[[TARGET]]

  n <- length(y)
  train_rows <- seq_len(floor(0.60 * n))
  val_rows <- seq(from = max(train_rows) + 1L, to = floor(0.80 * n))
  test_rows <- seq(from = max(val_rows) + 1L, to = n)
  X <- standardize_by_train(X, train_rows)

  varied <- apply(X[train_rows, , drop = FALSE], 2, stats::sd) > 1e-12
  X <- X[, varied, drop = FALSE]

  screen <- residual_screen(X, y, train_rows, val_rows, feature_set)
  fit_rows <- c(train_rows, val_rows)
  n_train <- length(train_rows)
  n_val <- length(val_rows)

  rows <- list()
  pred_rows <- list()
  for (method in METHODS) {
    cat(sprintf("Fitting %-7s | %s\n", feature_set, method))
    flush.console()
    fit <- gmdh_pmm(X[fit_rows, , drop = FALSE], y[fit_rows],
                    ctrl_for(method, n_train, n_val))
    pred <- predict(fit, X[test_rows, , drop = FALSE])

    rows[[length(rows) + 1L]] <- cbind(
      data.frame(
        target = TARGET,
        feature_set = feature_set,
        n_features = ncol(X),
        method = method,
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
      metrics(y[test_rows], pred),
      tail_metrics(y[test_rows], pred, 0.90, "top10_"),
      tail_metrics(y[test_rows], pred, 0.95, "top5_")
    )

    pred_rows[[length(pred_rows) + 1L]] <- data.frame(
      target = TARGET,
      feature_set = feature_set,
      method = method,
      time = d$time[test_rows],
      observed = y[test_rows],
      predicted = pred,
      residual = pred - y[test_rows],
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = do.call(rbind, rows),
    predictions = do.call(rbind, pred_rows),
    residual_screen = screen
  )
}

cat(sprintf("=== SRU time-aware C5 (%s, mode=%s, B=%d) ===\n\n", TARGET, FEATURE_MODE, B))
t0 <- Sys.time()
root <- find_root()
feature_sets <- if (FEATURE_MODE == "both") c("static", "dynamic") else FEATURE_MODE
parts <- lapply(feature_sets, function(fs) run_feature_set(root, fs))
summary <- do.call(rbind, lapply(parts, `[[`, "summary"))
summary <- summary[order(summary$feature_set, summary$nrmse), ]
predictions <- do.call(rbind, lapply(parts, `[[`, "predictions"))
screen <- do.call(rbind, lapply(parts, `[[`, "residual_screen"))

print(summary[, c("feature_set", "method", "nrmse", "mae", "bias",
                  "reserve_error", "top10_reserve_error", "top5_reserve_error",
                  "layers", "selected_nodes", "pmm2_share")],
      row.names = FALSE)
cat("\nResidual screen:\n")
print(screen[, c("feature_set", "dispatch", "gamma3", "gamma4", "g2", "val_nrmse")],
      row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
stem <- paste0("sru_timeaware_", TARGET)
utils::write.csv(summary, file.path(outdir, paste0(stem, "_summary.csv")), row.names = FALSE)
utils::write.csv(predictions, file.path(outdir, paste0(stem, "_predictions.csv")), row.names = FALSE)
utils::write.csv(screen, file.path(outdir, paste0(stem, "_residual_screen.csv")), row.names = FALSE)

png(file.path(outdir, paste0(stem, ".png")), width = 1120, height = 700, res = 115)
op <- par(mfrow = c(1, length(feature_sets)), mar = c(7, 4, 3, 1))
for (fs in feature_sets) {
  sub <- summary[summary$feature_set == fs, ]
  cols <- ifelse(sub$method == "PMM", "firebrick",
                 ifelse(sub$method == "LSE", "grey45",
                        ifelse(sub$method == "ridge-LSE", "steelblue",
                               ifelse(sub$method == "Huber", "darkgreen", "purple4"))))
  barplot(sub$nrmse, names.arg = sub$method, las = 2, col = cols,
          main = paste(TARGET, fs), ylab = "test NRMSE", cex.names = 0.75)
}
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s_{summary,predictions,residual_screen}.csv and %s.png | elapsed %.1fs\n",
            stem, stem, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
