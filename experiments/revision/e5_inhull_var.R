#!/usr/bin/env Rscript
# E5 -- IN-HULL VAR VARIANT (reviewer Q8: separate input-shift from coef-break).
#
# The original generate_var_break() in experiments/run_synthetic_var_break.R
# COUPLES the structural break to a post-break INPUT MEAN SHIFT
# (shift = c(2.0,-1.5,1.0)), pushing ~66% of the reserve rows OUTSIDE the
# training input hull. This confounds two distinct failure mechanisms:
#   (1) the y|x coefficient break (beta_post = beta_pre*(1+break_mag)), and
#   (2) the input-distribution shift (reserve rows extrapolate the train hull).
#
# This variant applies the SAME coefficient break but sets the post-break input
# mean shift to ZERO, so the VAR is input-stationary and reserve rows stay
# IN the training hull. Everything else identical: gamma3=1.5, break_mag=0.30,
# N=5000, gamma innovations, ReserveTW4 audit (reserve-aware, tail_weight=4).
#
# Hypothesis: with in-hull inputs, frac_outside_hull ~ 0 (vs ~0.66 original)
# and the catastrophic LSE upper tail (orig ReserveTW4 nrmse_max ~15.7)
# largely disappears -- the failure is driven by the input shift, not the break.
#
# Usage: Rscript experiments/revision/e5_inhull_var.R [N_SEEDS]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 30L
SEED0 <- 75001L
N <- 5000L
L_MAX <- 3L; F_KEEP <- 6L; TAIL_PROB <- 0.90; TAIL_WEIGHT <- 4; LEVEL <- 0.95
METHODS <- c("LSE", "PMM", "Huber")   # reviewer Q8 asks for these three

# ---- centred unit-variance gamma innovations (copied from the original) ----
rinnov <- function(n, gamma3 = 0) {
  if (abs(gamma3) < 1e-6) return(stats::rnorm(n))
  k <- 4 / gamma3^2
  x <- stats::rgamma(n, shape = k, rate = 1)
  z <- (x - k) / sqrt(k)
  if (gamma3 < 0) z <- -z
  z
}

# ---- IN-HULL generator: coefficient break kept, input mean shift = 0 ----
# Differs from generate_var_break() ONLY in `shift <- rep(0, d)`: the VAR(2)
# is then input-stationary, so the post-break / reserve inputs occupy the same
# region as the training inputs. The y|x coefficient jump is untouched.
generate_var_inhull <- function(seed, gamma3 = 1.5, break_mag = 0.30, N = 5000L) {
  set.seed(seed)
  d <- 3L
  A1 <- matrix(c(0.4, 0.1, 0.0, 0.0, 0.3, 0.1, 0.1, 0.0, 0.35), d, byrow = TRUE)
  A2 <- 0.2 * diag(d)
  X <- matrix(0, N, d)
  shift <- rep(0, d)                      # <<< IN-HULL: no post-break input shift
  for (t in 3:N) {
    mu <- if (t > N/2) shift else rep(0, d)
    X[t, ] <- mu + A1 %*% (X[t-1, ] - if (t-1 > N/2) shift else 0) +
                   A2 %*% (X[t-2, ] - if (t-2 > N/2) shift else 0) + stats::rnorm(d, sd = 0.5)
  }
  # y is an order-2 polynomial of x1,x2; coefficients jump at the break (KEPT).
  beta_pre  <- c(1.0, 0.8, -0.5, 0.6, -0.3, 0.2)   # 1, x1, x2, x1x2, x1^2, x2^2
  beta_post <- beta_pre * (1 + break_mag)
  feat <- function(x1, x2) cbind(1, x1, x2, x1*x2, x1^2, x2^2)
  pre <- seq_len(N) <= N/2
  y <- numeric(N)
  Fpre <- feat(X[,1], X[,2])
  y[pre]  <- Fpre[pre, ]  %*% beta_pre
  y[!pre] <- Fpre[!pre, ] %*% beta_post
  y <- y + rinnov(N, gamma3 = gamma3)
  list(X = X, y = y, transition = as.integer(N/2), n = N, gamma3 = gamma3,
       break_mag = break_mag)
}

# ---- helpers (copied from the original) ----
std <- function(X, ref) { mu <- colMeans(X[ref,,drop=FALSE]); sg <- apply(X[ref,,drop=FALSE],2,stats::sd)
  sg[!is.finite(sg)|sg<=1e-12] <- 1; sweep(sweep(X,2,mu,"-"),2,sg,"/") }
nrmse <- function(y,p){s<-stats::sd(y);if(!is.finite(s)||s<=0)s<-1;sqrt(mean((p-y)^2))/s}

ctrl_for <- function(m, ctype, n_tr, n_val, seed) {
  base <- list(L_max=L_MAX, F=F_KEEP, epsilon=-Inf, train_index=seq_len(n_tr),
               val_index=n_tr+seq_len(n_val), criterion=ctype, reserve_weight=1,
               tail_weight=TAIL_WEIGHT, tail_prob=TAIL_PROB, max_iter=60L)
  switch(m,
    PMM=do.call(gmdh_pmm_control,c(base,list(B=200L,force_method="auto",boot_seed=seed))),
    LSE=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="LSE"))),
    Huber=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="Huber"))))
}

ivh_diagnostics <- function(Xtr, Xte) {
  H <- Xtr %*% solve(crossprod(Xtr) + 1e-8*diag(ncol(Xtr))) %*% t(Xtr)
  lev_tr <- diag(H)
  mu <- colMeans(Xtr); S <- stats::cov(Xtr); Si <- solve(S + 1e-8*diag(ncol(Xtr)))
  maha <- function(M) mapply(function(i) { v <- M[i,]-mu; as.numeric(t(v)%*%Si%*%v) }, seq_len(nrow(M)))
  m_tr <- maha(Xtr); m_te <- maha(Xte)
  rng <- apply(Xtr, 2, range)
  outside <- mean(apply(Xte, 1, function(r) any(r < rng[1,] | r > rng[2,])))
  list(mean_lev_tr = mean(lev_tr), maha_tr_med = stats::median(m_tr),
       maha_te_med = stats::median(m_te), maha_ratio = stats::median(m_te)/stats::median(m_tr),
       frac_outside_hull = outside)
}

mshare <- function(fit, m) {
  tabs <- lapply(fit$layers, function(z) z$methods)
  num <- sum(vapply(tabs, function(t){v<-unname(t[m]); if(!length(v)||is.na(v))0 else v}, numeric(1)))
  den <- sum(vapply(tabs, sum, numeric(1))); if (den<=0) NA_real_ else num/den
}

# ---- ONE seed: ReserveTW4 audit only (criterion = reserve-aware) ----
run_audit_one <- function(seed, gamma3 = 1.5, break_mag = 0.30) {
  D <- generate_var_inhull(seed, gamma3, break_mag)
  T <- D$transition; n_tr <- floor(.6*T); n_val <- floor(.2*T)
  fr <- seq_len(n_tr+n_val); cal <- (n_tr+n_val+1L):T; test <- (T+1L):D$n
  Xs <- std(D$X, fr)
  ivh <- ivh_diagnostics(cbind(1,Xs[fr,,drop=FALSE]), cbind(1,Xs[test,,drop=FALSE]))
  rows <- list()
  ctype <- "reserve-aware"; clab <- paste0("ReserveTW",TAIL_WEIGHT)
  for (m in METHODS) {
    fit <- gmdh_pmm(Xs[fr,,drop=FALSE], D$y[fr], ctrl_for(m, ctype, n_tr, n_val, seed))
    pcal <- predict(fit, Xs[cal,,drop=FALSE]); pte <- predict(fit, Xs[test,,drop=FALSE])
    q <- as.numeric(stats::quantile(abs(D$y[cal]-pcal), LEVEL, names=FALSE))
    cov <- mean(D$y[test] >= pte-q & D$y[test] <= pte+q)
    rows[[length(rows)+1L]] <- data.frame(seed=seed, gamma3=gamma3, break_mag=break_mag,
      criterion=clab, method=m, nrmse_test=nrmse(D$y[test],pte),
      nrmse_cal=nrmse(D$y[cal],pcal), coverage=cov, width=2*q,
      pmm2_share=mshare(fit,"PMM2"), pmm3_share=mshare(fit,"PMM3"),
      frac_outside_hull=ivh$frac_outside_hull, maha_ratio=ivh$maha_ratio,
      stringsAsFactors=FALSE)
  }
  do.call(rbind, rows)
}

outdir <- "experiments/revision"; if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
t0 <- Sys.time()
cat(sprintf("=== E5 IN-HULL VAR audit (ReserveTW4) | seeds=%d gamma3=1.5 break=0.30 shift=0 ===\n", N_SEEDS))
all <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_audit_one(SEED0+s-1L)))
utils::write.csv(all, file.path(outdir,"e5_inhull_var_raw.csv"), row.names=FALSE)

agg <- do.call(rbind, by(all, list(all$method), function(g)
  data.frame(criterion=g$criterion[1], method=g$method[1], n_seed=nrow(g),
    nrmse_med=median(g$nrmse_test), nrmse_q1=quantile(g$nrmse_test,.25,names=FALSE),
    nrmse_q3=quantile(g$nrmse_test,.75,names=FALSE), nrmse_max=max(g$nrmse_test),
    cov_med=median(g$coverage), width_med=median(g$width),
    frac_outside_hull=median(g$frac_outside_hull), maha_ratio=median(g$maha_ratio),
    stringsAsFactors=FALSE)))
agg <- agg[order(agg$nrmse_med), ]
utils::write.csv(agg, file.path(outdir,"e5_inhull_var.csv"), row.names=FALSE)

cat(sprintf("\nReserve rows outside training hull (median frac_outside_hull) = %.4f\n",
            median(all$frac_outside_hull)))
cat(sprintf("(original coupled-shift generator: ~0.6594)\n\n"))
cat("=== E5 IN-HULL Reserve-test (median across seeds), ReserveTW4 ===\n")
print(agg[,c("method","n_seed","nrmse_med","nrmse_q1","nrmse_q3","nrmse_max",
             "cov_med","frac_outside_hull","maha_ratio")], row.names=FALSE, digits=4)

cat("\n=== ORIGINAL coupled-shift ReserveTW4 (for comparison) ===\n")
orig_path <- "experiments/results/synthetic_var_audit_summary.csv"
if (file.exists(orig_path)) {
  orig <- utils::read.csv(orig_path, stringsAsFactors=FALSE)
  orig <- orig[orig$criterion=="ReserveTW4" & orig$method %in% METHODS, ]
  orig <- orig[order(orig$nrmse_med), ]
  print(orig[,c("method","nrmse_med","nrmse_max","cov_med","frac_outside_hull","maha_ratio")],
        row.names=FALSE, digits=4)
} else cat("(original summary not found at", orig_path, ")\n")

cat(sprintf("\nSaved experiments/revision/e5_inhull_var{,_raw}.csv | %.0fs\n",
            as.numeric(difftime(Sys.time(),t0,"secs"))))
