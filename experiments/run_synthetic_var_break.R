#!/usr/bin/env Rscript
# Synthetic VAR(2) with a structural break -- controlled mechanism study for the
# ReserveTW4 LSE-GMDH failure (Paper 1 revision, Phase 3 / concern M3).
#
# Everything is known exactly here, so we can map the failure to its drivers:
#   (a) input-convex-hull (IVH) violation by the reserve-test rows,
#   (b) leverage / Cook's distance / Mahalanobis distance vs the training hull,
#   (c) the per-layer parameter / sample ratio (double-descent proxy).
#
# Innovations: we use centred, unit-variance GAMMA innovations rather than
# skew-Normal, because the skew-Normal skewness saturates at ~0.995 while the
# sweep needs gamma3 up to 2.0. Gamma(k) has skewness 2/sqrt(k), so a target
# gamma3 maps to k = 4 / gamma3^2 (gamma3 = 0 -> Gaussian).
#
# Usage: Rscript experiments/run_synthetic_var_break.R [N_SEEDS] [MODE]
#   MODE = "audit" (default, fixed gamma3=1.5/break=0.30) or "sweep".

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 30L
MODE    <- if (length(args) >= 2) args[2] else "audit"
SEED0 <- 75001L
N <- 5000L
L_MAX <- 3L; F_KEEP <- 6L; TAIL_PROB <- 0.90; TAIL_WEIGHT <- 4; LEVEL <- 0.95
METHODS <- c("PMM", "LSE", "ridge-LSE", "Huber", "L1")

# Centred unit-variance innovations.
#  innov = "gamma": asymmetric, target standardized skewness gamma3 (PMM2 regime).
#  innov = "platy": SYMMETRIC platykurtic via Beta(a,a), target excess kurtosis
#    gamma4 < 0 (PMM3 regime). Beta(a,a) has excess kurtosis -6/(2a+3), so
#    gamma4 -> a = -(3 + 6/gamma4)/2; uniform (a=1) gives gamma4 = -1.2.
rinnov <- function(n, gamma3 = 0, gamma4 = 0, innov = "gamma") {
  if (innov == "platy") {
    if (gamma4 >= -1e-6) return(stats::rnorm(n))
    a <- -(3 + 6 / gamma4) / 2
    a <- max(a, 0.05)
    x <- stats::rbeta(n, a, a)            # symmetric on (0,1), platykurtic
    (x - 0.5) / stats::sd(x)              # centre + unit variance
  } else {
    if (abs(gamma3) < 1e-6) return(stats::rnorm(n))
    k <- 4 / gamma3^2                     # gamma shape
    x <- stats::rgamma(n, shape = k, rate = 1)
    z <- (x - k) / sqrt(k)
    if (gamma3 < 0) z <- -z
    z
  }
}

# VAR(2) on d=3 latent predictors with a structural break at t = N/2:
# post-break the VAR mean shifts (pushing reserve rows outside the train hull)
# and the y|x coefficient map changes by +/- `break_mag`.
generate_var_break <- function(seed, gamma3 = 1.5, break_mag = 0.30, N = 5000L,
                               innov = "gamma", gamma4 = -1.2) {
  set.seed(seed)
  d <- 3L
  A1 <- matrix(c(0.4, 0.1, 0.0, 0.0, 0.3, 0.1, 0.1, 0.0, 0.35), d, byrow = TRUE)
  A2 <- 0.2 * diag(d)
  X <- matrix(0, N, d)
  shift <- c(2.0, -1.5, 1.0)              # post-break mean shift of the inputs
  for (t in 3:N) {
    mu <- if (t > N/2) shift else rep(0, d)
    X[t, ] <- mu + A1 %*% (X[t-1, ] - if (t-1 > N/2) shift else 0) +
                   A2 %*% (X[t-2, ] - if (t-2 > N/2) shift else 0) + stats::rnorm(d, sd = 0.5)
  }
  # y is an order-2 polynomial of x1,x2; coefficients jump at the break.
  beta_pre  <- c(1.0, 0.8, -0.5, 0.6, -0.3, 0.2)   # 1, x1, x2, x1x2, x1^2, x2^2
  beta_post <- beta_pre * (1 + break_mag)
  feat <- function(x1, x2) cbind(1, x1, x2, x1*x2, x1^2, x2^2)
  pre <- seq_len(N) <= N/2
  B <- ifelse(pre, 1, 0)
  y <- numeric(N)
  Fpre <- feat(X[,1], X[,2])
  y[pre]  <- Fpre[pre, ]  %*% beta_pre
  y[!pre] <- Fpre[!pre, ] %*% beta_post
  y <- y + rinnov(N, gamma3 = gamma3, gamma4 = gamma4, innov = innov)
  list(X = X, y = y, transition = as.integer(N/2), n = N, gamma3 = gamma3,
       break_mag = break_mag, innov = innov, gamma4 = gamma4)
}

std <- function(X, ref) { mu <- colMeans(X[ref,,drop=FALSE]); sg <- apply(X[ref,,drop=FALSE],2,stats::sd)
  sg[!is.finite(sg)|sg<=1e-12] <- 1; sweep(sweep(X,2,mu,"-"),2,sg,"/") }
nrmse <- function(y,p){s<-stats::sd(y);if(!is.finite(s)||s<=0)s<-1;sqrt(mean((p-y)^2))/s}

maxcond <- function(fit,Xf){out<-vector("list",length(fit$nodes));tr<-seq_len(fit$n_train);b<--Inf
  for(nd in fit$nodes){if(is.null(nd))next
    if(nd$type=="input"){out[[nd$id]]<-Xf[,nd$var];next}
    p1<-out[[nd$parents[1]]];p2<-out[[nd$parents[2]]];out[[nd$id]]<-kg2_predict(nd$theta,p1,p2)
    Z<-cbind(1,p1[tr],p2[tr],p1[tr]*p2[tr],p1[tr]^2,p2[tr]^2);sv<-svd(Z,nu=0,nv=0)$d
    cond<-if(length(sv)&&min(sv)>0)max(sv)/min(sv) else Inf;b<-max(b,log10(cond))};b}

ctrl_for <- function(m, ctype, n_tr, n_val, seed) {
  base <- list(L_max=L_MAX, F=F_KEEP, epsilon=-Inf, train_index=seq_len(n_tr),
               val_index=n_tr+seq_len(n_val), criterion=ctype, reserve_weight=1,
               tail_weight=TAIL_WEIGHT, tail_prob=TAIL_PROB, max_iter=60L)
  switch(m,
    PMM=do.call(gmdh_pmm_control,c(base,list(B=200L,force_method="auto",boot_seed=seed))),
    LSE=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="LSE"))),
    `ridge-LSE`=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="ridge-LSE",ridge_lambda=1e-2))),
    Huber=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="Huber"))),
    L1=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="L1"))))
}

# IVH / leverage diagnostics for the reserve-test rows vs the train hull.
ivh_diagnostics <- function(Xtr, Xte) {
  H <- Xtr %*% solve(crossprod(Xtr) + 1e-8*diag(ncol(Xtr))) %*% t(Xtr)
  lev_tr <- diag(H)
  mu <- colMeans(Xtr); S <- stats::cov(Xtr); Si <- solve(S + 1e-8*diag(ncol(Xtr)))
  maha <- function(M) mapply(function(i) { v <- M[i,]-mu; as.numeric(t(v)%*%Si%*%v) }, seq_len(nrow(M)))
  m_tr <- maha(Xtr); m_te <- maha(Xte)
  # fraction of reserve rows beyond the per-feature training range (hull proxy)
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

run_audit_one <- function(seed, gamma3 = 1.5, break_mag = 0.30,
                          innov = "gamma", gamma4 = -1.2) {
  D <- generate_var_break(seed, gamma3, break_mag, innov = innov, gamma4 = gamma4)
  T <- D$transition; n_tr <- floor(.6*T); n_val <- floor(.2*T)
  fr <- seq_len(n_tr+n_val); cal <- (n_tr+n_val+1L):T; test <- (T+1L):D$n
  Xs <- std(D$X, fr)
  ivh <- ivh_diagnostics(cbind(1,Xs[fr,,drop=FALSE]), cbind(1,Xs[test,,drop=FALSE]))
  rows <- list()
  for (ctype in c("MSE","reserve-aware")) {
    clab <- if (ctype=="MSE") "MSE" else paste0("ReserveTW",TAIL_WEIGHT)
    for (m in METHODS) {
      fit <- gmdh_pmm(Xs[fr,,drop=FALSE], D$y[fr], ctrl_for(m, ctype, n_tr, n_val, seed))
      pcal <- predict(fit, Xs[cal,,drop=FALSE]); pte <- predict(fit, Xs[test,,drop=FALSE])
      q <- as.numeric(stats::quantile(abs(D$y[cal]-pcal), LEVEL, names=FALSE))
      cov <- mean(D$y[test] >= pte-q & D$y[test] <= pte+q)
      rows[[length(rows)+1L]] <- data.frame(seed=seed, gamma3=gamma3, break_mag=break_mag,
        innov=innov, criterion=clab, method=m, nrmse_test=nrmse(D$y[test],pte),
        nrmse_cal=nrmse(D$y[cal],pcal), coverage=cov, width=2*q,
        max_log10_cond=maxcond(fit,Xs[fr,,drop=FALSE]),
        pmm2_share=mshare(fit,"PMM2"), pmm3_share=mshare(fit,"PMM3"),
        n_params=6*F_KEEP*fit$best$layer, n_train=n_tr,
        frac_outside_hull=ivh$frac_outside_hull, maha_ratio=ivh$maha_ratio,
        stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

outdir <- "experiments/results"; if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
t0 <- Sys.time()

if (MODE == "sweep") {
  cat("=== Synthetic VAR mechanism SWEEP ===\n")
  grid <- expand.grid(gamma3 = c(0,0.5,1.0,1.5,2.0), break_mag = c(0.10,0.20,0.30,0.50,1.00))
  all <- list()
  for (i in seq_len(nrow(grid))) {
    g <- grid$gamma3[i]; b <- grid$break_mag[i]
    for (s in seq_len(min(N_SEEDS,10L))) all[[length(all)+1L]] <- run_audit_one(SEED0+s-1L, g, b)
    cat(sprintf("  gamma3=%.1f break=%.2f done (%.0fs)\n", g, b, as.numeric(difftime(Sys.time(),t0,"secs"))))
  }
  raw <- do.call(rbind, all)
  utils::write.csv(raw, file.path(outdir,"synthetic_var_sweep_raw.csv"), row.names=FALSE)
  agg <- do.call(rbind, by(raw, list(raw$gamma3,raw$break_mag,raw$method,raw$criterion), function(g)
    data.frame(gamma3=g$gamma3[1],break_mag=g$break_mag[1],method=g$method[1],criterion=g$criterion[1],
      nrmse_med=median(g$nrmse_test),cov_med=median(g$coverage),
      frac_outside_hull=median(g$frac_outside_hull),maha_ratio=median(g$maha_ratio),stringsAsFactors=FALSE)))
  utils::write.csv(agg, file.path(outdir,"synthetic_var_sweep_summary.csv"), row.names=FALSE)
  cat("Saved synthetic_var_sweep_{raw,summary}.csv\n")
} else if (MODE == "platy") {
  # Symmetric platykurtic innovations (gamma4 = -1.2, gamma3 = 0) -> PMM3 regime.
  cat(sprintf("=== Synthetic VAR audit (PMM3 regime: symmetric platykurtic, gamma4=-1.2) | seeds=%d ===\n", N_SEEDS))
  all <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s)
    run_audit_one(SEED0+s-1L, gamma3 = 0, break_mag = 0.30, innov = "platy", gamma4 = -1.2)))
  utils::write.csv(all, file.path(outdir,"synthetic_var_platy_raw.csv"), row.names=FALSE)
  agg <- do.call(rbind, by(all, list(all$criterion,all$method), function(g)
    data.frame(criterion=g$criterion[1],method=g$method[1],n_seed=nrow(g),
      nrmse_med=median(g$nrmse_test),nrmse_q1=quantile(g$nrmse_test,.25,names=FALSE),
      nrmse_q3=quantile(g$nrmse_test,.75,names=FALSE),nrmse_max=max(g$nrmse_test),
      cov_med=median(g$coverage),pmm2_med=median(g$pmm2_share,na.rm=TRUE),
      pmm3_med=median(g$pmm3_share,na.rm=TRUE),cond_med=median(g$max_log10_cond),stringsAsFactors=FALSE)))
  agg <- agg[order(agg$criterion, agg$nrmse_med), ]
  utils::write.csv(agg, file.path(outdir,"synthetic_var_platy_summary.csv"), row.names=FALSE)
  cat("\n=== Reserve-test (median across seeds), PMM3 regime ===\n")
  print(agg[,c("criterion","method","n_seed","nrmse_med","nrmse_q1","nrmse_q3","nrmse_max",
               "cov_med","pmm2_med","pmm3_med")], row.names=FALSE, digits=4)
  cat(sprintf("\nSaved synthetic_var_platy_{raw,summary}.csv | %.0fs\n", as.numeric(difftime(Sys.time(),t0,"secs"))))
} else {
  cat(sprintf("=== Synthetic VAR audit | seeds=%d gamma3=1.5 break=0.30 ===\n", N_SEEDS))
  all <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_audit_one(SEED0+s-1L)))
  utils::write.csv(all, file.path(outdir,"synthetic_var_audit_raw.csv"), row.names=FALSE)
  agg <- do.call(rbind, by(all, list(all$criterion,all$method), function(g)
    data.frame(criterion=g$criterion[1],method=g$method[1],n_seed=nrow(g),
      nrmse_med=median(g$nrmse_test),nrmse_q1=quantile(g$nrmse_test,.25,names=FALSE),
      nrmse_q3=quantile(g$nrmse_test,.75,names=FALSE),nrmse_max=max(g$nrmse_test),
      cov_med=median(g$coverage),width_med=median(g$width),cond_med=median(g$max_log10_cond),
      frac_outside_hull=median(g$frac_outside_hull),maha_ratio=median(g$maha_ratio),stringsAsFactors=FALSE)))
  agg <- agg[order(agg$criterion, agg$nrmse_med), ]
  utils::write.csv(agg, file.path(outdir,"synthetic_var_audit_summary.csv"), row.names=FALSE)
  cat("\n=== Reserve-test (median across seeds) ===\n")
  print(agg[,c("criterion","method","n_seed","nrmse_med","nrmse_q1","nrmse_q3","nrmse_max",
               "cov_med","cond_med","frac_outside_hull","maha_ratio")], row.names=FALSE, digits=4)
  cat(sprintf("\nSaved synthetic_var_audit_{raw,summary}.csv | %.0fs\n", as.numeric(difftime(Sys.time(),t0,"secs"))))
}
