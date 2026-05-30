#!/usr/bin/env Rscript
# E8 -- ReserveTW4 selection-criterion sensitivity (Paper 1 revision).
#
# Question: is the qualitative ranking under the reserve-aware selection
# criterion (PMM best median, LSE catastrophic upper tail) ROBUST to the two
# free knobs of the criterion -- the tail weight w_T (tail_weight) and the
# decile cutoff q (tail_prob)?
#
# We re-run the 4-split reserve audit of run_synthetic_var_break.R at the
# headline config (gamma3=1.5, break_mag=0.30, gamma innovations), over a modest
# 15 seeds, for the cross of
#     tail_weight in {2, 4, 8}  x  tail_prob in {0.85, 0.90, 0.95}   (9 cells)
# restricting the inner estimators to {PMM, LSE, Huber} and the selection
# criterion to "reserve-aware". For each cell we report, per method, the median
# reserve-test NRMSE; the LSE worst-case (max) NRMSE; and two booleans:
#   pmm_lowest_median -- is PMM the lowest-median method in the cell?
#   lse_largest_max   -- does LSE carry the largest max NRMSE in the cell?
#
# HEADLINE TO TEST: the ranking is INSENSITIVE to w_T and q (PMM best in all 9).
#
# Usage: Rscript experiments/revision/e8_reservetw4_sensitivity.R [N_SEEDS]
#   default N_SEEDS = 15.  Audit at N=5000 costs ~1-2 s / seed / cell.

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args    <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 15L
SEED0   <- 75001L
N       <- 5000L
L_MAX   <- 3L
F_KEEP  <- 6L
LEVEL   <- 0.95
GAMMA3  <- 1.5
BREAK   <- 0.30
INNOV   <- "gamma"
GAMMA4  <- -1.2           # unused under gamma innovations; kept for the generator signature
METHODS <- c("PMM", "LSE", "Huber")

# Sensitivity grid for the ReserveTW selection criterion.
TW_GRID <- c(2, 4, 8)            # tail_weight  (w_T): default 4
TP_GRID <- c(0.85, 0.90, 0.95)   # tail_prob    (q):   default 0.90

# ---- innovations + generator (verbatim from run_synthetic_var_break.R) -------
rinnov <- function(n, gamma3 = 0, gamma4 = 0, innov = "gamma") {
  if (innov == "platy") {
    if (gamma4 >= -1e-6) return(stats::rnorm(n))
    a <- -(3 + 6 / gamma4) / 2
    a <- max(a, 0.05)
    x <- stats::rbeta(n, a, a)
    (x - 0.5) / stats::sd(x)
  } else {
    if (abs(gamma3) < 1e-6) return(stats::rnorm(n))
    k <- 4 / gamma3^2
    x <- stats::rgamma(n, shape = k, rate = 1)
    z <- (x - k) / sqrt(k)
    if (gamma3 < 0) z <- -z
    z
  }
}

generate_var_break <- function(seed, gamma3 = 1.5, break_mag = 0.30, N = 5000L,
                               innov = "gamma", gamma4 = -1.2) {
  set.seed(seed)
  d <- 3L
  A1 <- matrix(c(0.4, 0.1, 0.0, 0.0, 0.3, 0.1, 0.1, 0.0, 0.35), d, byrow = TRUE)
  A2 <- 0.2 * diag(d)
  X <- matrix(0, N, d)
  shift <- c(2.0, -1.5, 1.0)
  for (t in 3:N) {
    mu <- if (t > N/2) shift else rep(0, d)
    X[t, ] <- mu + A1 %*% (X[t-1, ] - if (t-1 > N/2) shift else 0) +
                   A2 %*% (X[t-2, ] - if (t-2 > N/2) shift else 0) + stats::rnorm(d, sd = 0.5)
  }
  beta_pre  <- c(1.0, 0.8, -0.5, 0.6, -0.3, 0.2)
  beta_post <- beta_pre * (1 + break_mag)
  feat <- function(x1, x2) cbind(1, x1, x2, x1*x2, x1^2, x2^2)
  pre <- seq_len(N) <= N/2
  y <- numeric(N)
  Fpre <- feat(X[,1], X[,2])
  y[pre]  <- Fpre[pre, ]  %*% beta_pre
  y[!pre] <- Fpre[!pre, ] %*% beta_post
  y <- y + rinnov(N, gamma3 = gamma3, gamma4 = gamma4, innov = innov)
  list(X = X, y = y, transition = as.integer(N/2), n = N, gamma3 = gamma3,
       break_mag = break_mag, innov = innov, gamma4 = gamma4)
}

std <- function(X, ref) {
  mu <- colMeans(X[ref,,drop=FALSE]); sg <- apply(X[ref,,drop=FALSE], 2, stats::sd)
  sg[!is.finite(sg) | sg <= 1e-12] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sg, "/")
}
nrmse <- function(y, p) { s <- stats::sd(y); if (!is.finite(s) || s <= 0) s <- 1; sqrt(mean((p-y)^2)) / s }

# Control builder parametrized by tail_weight / tail_prob (reserve-aware only).
ctrl_for <- function(m, n_tr, n_val, seed, tw, tp) {
  base <- list(L_max = L_MAX, F = F_KEEP, epsilon = -Inf,
               train_index = seq_len(n_tr), val_index = n_tr + seq_len(n_val),
               criterion = "reserve-aware", reserve_weight = 1,
               tail_weight = tw, tail_prob = tp, max_iter = 60L)
  switch(m,
    PMM   = do.call(gmdh_pmm_control, c(base, list(B = 200L, force_method = "auto", boot_seed = seed))),
    LSE   = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "LSE"))),
    Huber = do.call(gmdh_pmm_control, c(base, list(B = 0, force_method = "Huber"))))
}

# One reserve audit (one seed, one cell): returns per-method reserve NRMSE.
run_audit_one <- function(seed, tw, tp) {
  D <- generate_var_break(seed, GAMMA3, BREAK, N = N, innov = INNOV, gamma4 = GAMMA4)
  T <- D$transition
  n_tr  <- floor(.6 * T); n_val <- floor(.2 * T)
  fr   <- seq_len(n_tr + n_val)            # pre-break train (60%) + selection-val (20%)
  cal  <- (n_tr + n_val + 1L):T            # interval-calibration (20% of pre-break)
  test <- (T + 1L):D$n                     # entire post-break block = reserve test
  Xs <- std(D$X, fr)
  out <- lapply(METHODS, function(m) {
    fit <- gmdh_pmm(Xs[fr,,drop=FALSE], D$y[fr], ctrl_for(m, n_tr, n_val, seed, tw, tp))
    pte <- predict(fit, Xs[test,,drop=FALSE])
    data.frame(seed = seed, tail_weight = tw, tail_prob = tp, method = m,
               nrmse_test = nrmse(D$y[test], pte), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# ---- run the 9-cell x 15-seed grid -------------------------------------------
outdir <- "experiments/revision"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
t0 <- Sys.time()
cat(sprintf("=== E8 ReserveTW sensitivity | seeds=%d  gamma3=%.1f break=%.2f innov=%s ===\n",
            N_SEEDS, GAMMA3, BREAK, INNOV))
cat(sprintf("    grid: tail_weight in {%s} x tail_prob in {%s}  (9 cells), methods {%s}\n",
            paste(TW_GRID, collapse=", "), paste(TP_GRID, collapse=", "),
            paste(METHODS, collapse=", ")))

raw <- list()
for (tw in TW_GRID) for (tp in TP_GRID) {
  for (s in seq_len(N_SEEDS)) raw[[length(raw)+1L]] <- run_audit_one(SEED0 + s - 1L, tw, tp)
  cat(sprintf("  tail_weight=%g tail_prob=%.2f done (%.0fs)\n",
              tw, tp, as.numeric(difftime(Sys.time(), t0, "secs"))))
}
raw <- do.call(rbind, raw)
utils::write.csv(raw, file.path(outdir, "e8_reservetw4_sensitivity_raw.csv"), row.names = FALSE)

# ---- per-cell summary --------------------------------------------------------
cells <- list()
for (tw in TW_GRID) for (tp in TP_GRID) {
  sub <- raw[raw$tail_weight == tw & raw$tail_prob == tp, ]
  med <- tapply(sub$nrmse_test, sub$method, stats::median)
  mx  <- tapply(sub$nrmse_test, sub$method, max)
  med <- med[METHODS]; mx <- mx[METHODS]               # fix method order
  lowest_method   <- names(med)[which.min(med)]
  largest_max_mtd <- names(mx)[which.max(mx)]
  cells[[length(cells)+1L]] <- data.frame(
    tail_weight = tw, tail_prob = tp, n_seed = N_SEEDS,
    pmm_med   = unname(med["PMM"]),
    lse_med   = unname(med["LSE"]),
    huber_med = unname(med["Huber"]),
    lse_max   = unname(mx["LSE"]),
    pmm_max   = unname(mx["PMM"]),
    huber_max = unname(mx["Huber"]),
    pmm_lowest_median = (lowest_method == "PMM"),
    lse_largest_max   = (largest_max_mtd == "LSE"),
    stringsAsFactors = FALSE)
}
cells <- do.call(rbind, cells)
utils::write.csv(cells, file.path(outdir, "e8_reservetw4_sensitivity.csv"), row.names = FALSE)

# ---- verdict -----------------------------------------------------------------
pmm_all  <- all(cells$pmm_lowest_median)
lse_all  <- all(cells$lse_largest_max)

cat("\n=== Per-cell summary (reserve-test NRMSE) ===\n")
print(cells[, c("tail_weight","tail_prob","pmm_med","lse_med","huber_med",
                "lse_max","pmm_lowest_median","lse_largest_max")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nPMM lowest-median in %d/%d cells.\n", sum(cells$pmm_lowest_median), nrow(cells)))
cat(sprintf("LSE largest-max  in %d/%d cells.\n", sum(cells$lse_largest_max), nrow(cells)))
cat(sprintf("\nVERDICT  PMM-best-in-all-cells: %s\n", if (pmm_all) "YES" else "NO"))
cat(sprintf("         LSE-largest-max-in-all-cells: %s\n", if (lse_all) "YES" else "NO"))
cat(sprintf("\nSaved e8_reservetw4_sensitivity{,_raw}.csv | %.0fs\n",
            as.numeric(difftime(Sys.time(), t0, "secs"))))
