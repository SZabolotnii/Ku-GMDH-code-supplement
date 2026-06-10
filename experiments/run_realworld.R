#!/usr/bin/env Rscript
# C5 real-world / built-in candidate comparison for Paper 1.
#
# This script takes the real-data shortlist and runs the full GMDH tournament
# with identical outer-loop settings for:
# PMM-auto, LSE, ridge-LSE, Huber, and L1.
#
# It is still a candidate-stage experiment, not a final flagship. The goal is
# to decide which screened datasets deserve deeper C5 treatment.
#
# Run from the package directory:
#   Rscript experiments/run_realworld.R [candidate|all|domain|extended|insurance|external] [R] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
CAND <- if (length(args) >= 1) args[1] else "airquality_ozone"
R_REP <- if (length(args) >= 2) as.integer(args[2]) else 50L
B <- if (length(args) >= 3) as.integer(args[3]) else 80L
TEST_FRAC <- 0.25
L_MAX <- 3L
F_KEEP <- 6L

load_data <- function(pkg, name) {
  env <- new.env(parent = globalenv())
  suppressWarnings(utils::data(list = name, package = pkg, envir = env))
  if (!exists(name, envir = env, inherits = FALSE)) {
    stop("Dataset ", pkg, "::", name, " is unavailable in this R installation.")
  }
  get(name, envir = env, inherits = FALSE)
}

model_x <- function(formula, data) {
  X <- stats::model.matrix(formula, data = data)
  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]
  finite <- apply(X, 2, function(z) all(is.finite(z)))
  varied <- apply(X, 2, function(z) stats::sd(z) > 1e-12)
  X <- X[, finite & varied, drop = FALSE]
  storage.mode(X) <- "double"
  X
}

find_root <- function() {
  here <- normalizePath(getwd(), mustWork = TRUE)
  cur <- here
  repeat {
    if (file.exists(file.path(cur, "DESCRIPTION")) &&
        dir.exists(file.path(cur, "R"))) return(cur)
    if (file.exists(file.path(cur, "ROADMAP.md")) &&
        file.exists(file.path(cur, "paper-1-gmdh-pmm", "code", "DESCRIPTION"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Cannot locate repository root from ", here)
    cur <- parent
  }
}

root <- find_root()
standalone_repo <- file.exists(file.path(root, "DESCRIPTION")) &&
  dir.exists(file.path(root, "R"))

manifest_candidates <- if (standalone_repo) {
  c(file.path(root, "datasets", "external", "external_candidates.csv"),
    file.path(root, "shared", "datasets", "external", "external_candidates.csv"))
} else {
  c(file.path(root, "shared", "datasets", "external", "external_candidates.csv"),
    file.path(root, "datasets", "external", "external_candidates.csv"))
}
external_manifest_path <- manifest_candidates[file.exists(manifest_candidates)][1]
if (!length(external_manifest_path) || is.na(external_manifest_path)) {
  external_manifest_path <- character(0)
}

external_manifest <- function() {
  if (!length(external_manifest_path) || !file.exists(external_manifest_path)) return(data.frame())
  utils::read.csv(external_manifest_path, stringsAsFactors = FALSE)
}

resolve_external_path <- function(path, must_work = TRUE) {
  if (grepl("^(/|[A-Za-z]:)", path)) {
    candidates <- path
  } else {
    candidates <- c(file.path(root, path))
    if (standalone_repo) {
      candidates <- c(
        candidates,
        file.path(root, "datasets", "external", basename(path)),
        file.path(root, sub("^shared/datasets/external/processed/", "datasets/external/", path))
      )
    }
  }
  candidates <- unique(candidates)
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) && !is.na(hit)) return(hit)
  if (must_work) {
    stop("External candidate file not found. Tried: ",
         paste(candidates, collapse = "; "),
         "\nOnly selected public datasets are shipped in this supplement; ",
         "see datasets/README.md.")
  }
  NULL
}

external_candidate_ids <- function() {
  x <- external_manifest()
  if (!nrow(x)) return(character(0))
  if (!standalone_repo) return(x$id)
  keep <- vapply(x$path, function(path) !is.null(resolve_external_path(path, must_work = FALSE)),
                 logical(1))
  x$id[keep]
}

external_candidate_data <- function(id) {
  manifest <- external_manifest()
  if (!nrow(manifest) || !(id %in% manifest$id)) return(NULL)
  row <- manifest[match(id, manifest$id), ]
  data_path <- resolve_external_path(row$path[[1]])
  d <- stats::na.omit(utils::read.csv(data_path, stringsAsFactors = TRUE))
  formula <- stats::as.formula(row$formula[[1]])
  mf <- stats::model.frame(formula, data = d, na.action = stats::na.omit)
  y <- stats::model.response(mf)
  list(
    id = id,
    label = row$label[[1]],
    y = y,
    X = model_x(formula, mf)
  )
}

candidate_data <- function(id) {
  ext <- external_candidate_data(id)
  if (!is.null(ext)) return(ext)

  data(airquality, package = "datasets")
  data(iris, package = "datasets")
  data(ToothGrowth, package = "datasets")

  if (id == "airquality_ozone") {
    d <- stats::na.omit(airquality[, c("Ozone", "Solar.R", "Wind", "Temp")])
    return(list(
      id = id,
      label = "Airquality: Ozone ~ Solar.R + Wind + Temp",
      y = d$Ozone,
      X = as.matrix(d[, c("Solar.R", "Wind", "Temp")])
    ))
  }
  if (id == "iris_versicolor") {
    d <- subset(iris, Species == "versicolor")
    return(list(
      id = id,
      label = "Iris versicolor: Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width",
      y = d$Sepal.Length,
      X = as.matrix(d[, c("Sepal.Width", "Petal.Length", "Petal.Width")])
    ))
  }
  if (id == "toothgrowth") {
    d <- ToothGrowth
    X <- cbind(dose = d$dose, suppVC = as.numeric(d$supp == "VC"))
    return(list(
      id = id,
      label = "ToothGrowth: len ~ dose + supp",
      y = d$len,
      X = as.matrix(X)
    ))
  }
  if (id == "boot_acme_returns") {
    acme <- load_data("boot", "acme")
    d <- data.frame(
      acme = acme$acme[-1],
      market = acme$market[-1],
      market_lag1 = acme$market[-nrow(acme)],
      acme_lag1 = acme$acme[-nrow(acme)]
    )
    return(list(
      id = id,
      label = "Acme monthly excess returns: acme_t ~ market_t + market_{t-1} + acme_{t-1}",
      y = d$acme,
      X = as.matrix(d[, c("market", "market_lag1", "acme_lag1")])
    ))
  }
  if (id == "islr_credit_balance") {
    d <- stats::na.omit(load_data("ISLR2", "Credit"))
    return(list(
      id = id,
      label = "Credit balance ~ account covariates",
      y = d$Balance,
      X = model_x(~ Income + Limit + Rating + Cards + Age + Education + Own + Student + Married + Region, d)
    ))
  }
  if (id == "survival_solder_skips") {
    d <- stats::na.omit(load_data("survival", "solder"))
    return(list(
      id = id,
      label = "Soldering experiment: skips ~ process factors",
      y = d$skips,
      X = model_x(~ Opening + Solder + Mask + PadType + Panel, d)
    ))
  }
  if (id == "boot_motor_accel") {
    d <- stats::na.omit(load_data("boot", "motor"))
    X <- cbind(times = d$times, times2 = d$times^2, times3 = d$times^3)
    return(list(
      id = id,
      label = "Motorcycle crash: accel ~ polynomial time basis",
      y = d$accel,
      X = X
    ))
  }
  if (id == "mass_cars93_mpg") {
    d <- stats::na.omit(load_data("MASS", "Cars93"))
    return(list(
      id = id,
      label = "Cars93: MPG.city ~ vehicle covariates",
      y = d$MPG.city,
      X = model_x(~ Horsepower + EngineSize + RPM + Weight + Type + Origin, d)
    ))
  }
  if (id == "insurance_autobi_loss") {
    d <- stats::na.omit(load_data("insuranceData", "AutoBi"))
    d$ATTORNEY <- factor(d$ATTORNEY)
    d$CLMSEX <- factor(d$CLMSEX)
    d$MARITAL <- factor(d$MARITAL)
    d$CLMINSUR <- factor(d$CLMINSUR)
    d$SEATBELT <- factor(d$SEATBELT)
    return(list(
      id = id,
      label = "AutoBi bodily injury: LOSS ~ claimant/legal covariates",
      y = d$LOSS,
      X = model_x(~ ATTORNEY + CLMSEX + MARITAL + CLMINSUR + SEATBELT + CLMAGE, d)
    ))
  }
  if (id == "insurance_autoclaims_paid") {
    d <- stats::na.omit(load_data("insuranceData", "AutoClaims"))
    return(list(
      id = id,
      label = "AutoClaims property claims: PAID ~ rating covariates",
      y = d$PAID,
      X = model_x(~ STATE + CLASS + GENDER + AGE, d)
    ))
  }
  if (id == "islr_wage") {
    d <- stats::na.omit(load_data("ISLR2", "Wage"))
    return(list(
      id = id,
      label = "Wage ~ demographic and job covariates",
      y = d$wage,
      X = model_x(~ year + age + maritl + race + education + jobclass + health + health_ins, d)
    ))
  }
  if (id == "mass_boston_medv") {
    d <- stats::na.omit(load_data("MASS", "Boston"))
    return(list(
      id = id,
      label = "Boston housing: medv ~ structural covariates",
      y = d$medv,
      X = model_x(~ crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + lstat, d)
    ))
  }
  stop("Unknown candidate: ", id)
}

ctrl_for <- function(method, seed) {
  base <- list(L_max = L_MAX, F = F_KEEP, epsilon = -Inf, seed = seed)
  switch(method,
         PMM = do.call(gmdh_pmm_control, c(base, list(B = B, force_method = "auto"))),
         LSE = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "LSE"))),
         `ridge-LSE` = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "ridge-LSE"))),
         Huber = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "Huber"))),
         L1 = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "L1"))),
         stop("Unknown method: ", method))
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

one_rep <- function(dat, rep_id, seed) {
  set.seed(seed)
  n <- length(dat$y)
  n_test <- max(10L, ceiling(TEST_FRAC * n))
  test <- sort(sample.int(n, n_test))
  train <- setdiff(seq_len(n), test)
  xs <- standardize_split(dat$X, train, test)
  y_train <- dat$y[train]
  y_test <- dat$y[test]
  y_scale <- stats::sd(y_test)
  if (!is.finite(y_scale) || y_scale <= 1e-12) y_scale <- stats::sd(y_train)
  if (!is.finite(y_scale) || y_scale <= 1e-12) y_scale <- 1

  rows <- list()
  for (method in c("PMM", "LSE", "ridge-LSE", "Huber", "L1")) {
    fit <- try(gmdh_pmm(xs$X_train, y_train, ctrl_for(method, seed + match(method, c("PMM", "LSE", "ridge-LSE", "Huber", "L1")))),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      rows[[length(rows) + 1]] <- data.frame(
        candidate = dat$id, rep = rep_id, method = method,
        nrmse = NA_real_, mae = NA_real_, bias = NA_real_,
        pmm2_share = NA_real_, pmm3_share = NA_real_, layers = NA_integer_,
        ok = FALSE, stringsAsFactors = FALSE)
      next
    }
    pred <- predict(fit, xs$X_test)
    e <- pred - y_test
    rows[[length(rows) + 1]] <- data.frame(
      candidate = dat$id,
      rep = rep_id,
      method = method,
      nrmse = sqrt(mean(e^2)) / y_scale,
      mae = mean(abs(e)),
      bias = mean(e),
      pmm2_share = method_share(fit, "PMM2"),
      pmm3_share = method_share(fit, "PMM3"),
      layers = length(fit$layers),
      ok = TRUE,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  best <- out$method[which.min(out$nrmse)]
  out$winner <- out$method == best
  out
}

summarize <- function(raw) {
  parts <- split(raw[raw$ok, ], list(raw$candidate[raw$ok], raw$method[raw$ok]), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(g) {
    data.frame(
      candidate = g$candidate[1],
      method = g$method[1],
      nrmse_mean = round(mean(g$nrmse), 3),
      nrmse_sd = round(stats::sd(g$nrmse), 3),
      mae_mean = round(mean(g$mae), 3),
      bias_mean = round(mean(g$bias), 3),
      win_rate = round(mean(g$winner), 3),
      pmm2_share = round(mean(g$pmm2_share, na.rm = TRUE), 3),
      pmm3_share = round(mean(g$pmm3_share, na.rm = TRUE), 3),
      layers_mean = round(mean(g$layers), 2),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$candidate, out$nrmse_mean), ]
}

builtin_candidates <- c("airquality_ozone", "iris_versicolor", "toothgrowth")
domain_candidates <- c("boot_acme_returns", "islr_credit_balance", "survival_solder_skips",
                       "boot_motor_accel", "mass_cars93_mpg")
insurance_candidates <- c("insurance_autobi_loss", "insurance_autoclaims_paid")
extended_candidates <- c(domain_candidates, "islr_wage", "mass_boston_medv")
external_candidates <- external_candidate_ids()
candidates <- switch(CAND,
                     all = builtin_candidates,
                     builtins = builtin_candidates,
                     domain = domain_candidates,
                     insurance = insurance_candidates,
                     extended = extended_candidates,
                     external = external_candidates,
                     CAND)
if (!length(candidates)) stop("No candidates resolved for mode: ", CAND)
cat(sprintf("=== C5 GMDH real-world candidate run (candidate=%s, R=%d, B=%d) ===\n\n",
            paste(candidates, collapse = ","), R_REP, B))
t0 <- Sys.time()
raw <- do.call(rbind, unlist(lapply(seq_along(candidates), function(i) {
  dat <- candidate_data(candidates[i])
  lapply(seq_len(R_REP), function(r) one_rep(dat, r, 30000L + i * 1000L + r))
}), recursive = FALSE))
res <- summarize(raw)
print(res, row.names = FALSE)

outdir <- "experiments/results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
suffix <- switch(CAND,
                 all = "",
                 builtins = "",
                 domain = "_domain",
                 insurance = "_insurance",
                 extended = "_extended",
                 external = "_external",
                 paste0("_", gsub("[^A-Za-z0-9]+", "_", CAND)))
stem <- paste0("realworld_c5", suffix)
utils::write.csv(raw, file.path(outdir, paste0(stem, "_raw.csv")), row.names = FALSE)
utils::write.csv(res, file.path(outdir, paste0(stem, "_summary.csv")), row.names = FALSE)

png(file.path(outdir, paste0(stem, ".png")),
    width = max(920, 280 * length(candidates)), height = 520, res = 110)
op <- par(mfrow = c(1, length(candidates)), mar = c(6, 4, 3, 1))
for (id in candidates) {
  sub <- res[res$candidate == id, ]
  cols <- ifelse(sub$method == "PMM", "firebrick",
                 ifelse(sub$method == "LSE", "grey45",
                        ifelse(sub$method == "ridge-LSE", "steelblue",
                               ifelse(sub$method == "Huber", "darkgreen", "purple4"))))
  barplot(sub$nrmse_mean, names.arg = sub$method, las = 2, col = cols,
          ylab = "test NRMSE", main = id, cex.names = 0.75)
}
par(op); invisible(dev.off())

cat(sprintf("\nSaved %s_summary.csv, %s_raw.csv, %s.png | elapsed %.1fs\n",
            stem, stem, stem,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
