#!/usr/bin/env Rscript
# Volve layer-by-layer survival + forced-survival-set control experiment
# (Paper 1 revision, Phase 1.7 / concern M3).
#
# Question M3: is the LSE-GMDH reserve-window failure SELECTION-induced (LSE
# selects a worse structure) or ESTIMATOR-induced (same structure, worse
# coefficients)? We (1) record the per-layer method shares and the best-node
# parent pair for each estimator, then (2) freeze the best parent pair chosen
# under one estimator and re-fit that SAME KG-2 partial model with every
# estimator, comparing reserve-test NRMSE. If the failure persists under a
# frozen structure, it is estimator-induced.

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

L_MAX <- 3L; F_KEEP <- 6L; TAIL_PROB <- 0.90; TAIL_WEIGHT <- 4
METHODS <- c("PMM", "LSE", "ridge-LSE", "Huber", "L1")
SEED <- 75001L

std <- function(X, ref){mu<-colMeans(X[ref,,drop=FALSE]);sg<-apply(X[ref,,drop=FALSE],2,stats::sd)
  sg[!is.finite(sg)|sg<=1e-12]<-1;sweep(sweep(X,2,mu,"-"),2,sg,"/")}
nrmse <- function(y,p){s<-stats::sd(y);if(!is.finite(s)||s<=0)s<-1;sqrt(mean((p-y)^2))/s}

ctrl_for <- function(m,n_tr,n_val){
  base<-list(L_max=L_MAX,F=F_KEEP,epsilon=-Inf,train_index=seq_len(n_tr),
    val_index=n_tr+seq_len(n_val),criterion="reserve-aware",reserve_weight=1,
    tail_weight=TAIL_WEIGHT,tail_prob=TAIL_PROB,max_iter=60L)
  switch(m,
    PMM=do.call(gmdh_pmm_control,c(base,list(B=200L,force_method="auto",boot_seed=SEED))),
    LSE=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="LSE"))),
    `ridge-LSE`=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="ridge-LSE",ridge_lambda=1e-2))),
    Huber=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="Huber"))),
    L1=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="L1"))))
}

D <- load_volve_drilling(segment="gr", target="log_rop5")
T <- D$transition; n_tr<-floor(.6*T); n_val<-floor(.2*T)
fr<-seq_len(n_tr+n_val); test<-(T+1L):D$n
Xs<-std(D$X,fr)

cat("=== Layer-by-layer survival & best-node by estimator (reserve-aware) ===\n")
best_pairs <- list()
for (m in METHODS) {
  fit <- gmdh_pmm(Xs[fr,,drop=FALSE], D$y[fr], ctrl_for(m,n_tr,n_val))
  bn <- fit$nodes[[fit$best_id]]
  shares <- sapply(fit$layers, function(z) paste(names(z$methods), z$methods, sep="=", collapse=","))
  pred <- predict(fit, Xs[test,,drop=FALSE])
  best_pairs[[m]] <- bn$parents
  bm <- if (is.null(bn$method)) "input" else bn$method
  cat(sprintf("%-10s best_layer=%d best_node parents=(%s) method=%s | reserve NRMSE=%.3f\n  layer methods: %s\n",
      m, fit$best$layer, paste(bn$parents,collapse=","), bm,
      nrmse(D$y[test],pred), paste(shares, collapse=" | ")))
}

# --- Forced-survival-set control: freeze best parent pair, swap estimator ---
cat("\n=== Control: freeze best pair, re-estimate with each inner method ===\n")
cat("Predictor index: 1=swob 2=rpm 3=tqa 4=dept (layer-1 pairs reference inputs)\n")
fit_pair <- function(a, b, method) {
  v1tr<-Xs[fr,a]; v2tr<-Xs[fr,b]; v1te<-Xs[test,a]; v2te<-Xs[test,b]
  ctrl<-ctrl_for(method,n_tr,n_val)
  est<-inner_estimate(v1tr[seq_len(n_tr)], v2tr[seq_len(n_tr)], D$y[fr][seq_len(n_tr)], ctrl)
  pred<-kg2_predict(est$theta, v1te, v2te)
  list(nrmse=nrmse(D$y[test],pred), method=est$method, theta=est$theta)
}
for (src in c("PMM","LSE")) {
  pr <- best_pairs[[src]]
  if (length(pr)!=2 || any(pr>4)) { cat(sprintf("(%s best node is multi-layer; skipping pair freeze)\n",src)); next }
  cat(sprintf("\n-- Frozen structure = best pair from %s: (%d,%d) --\n", src, pr[1], pr[2]))
  for (m in c("LSE","PMM","Huber","L1")) {
    r <- fit_pair(pr[1], pr[2], m)
    cat(sprintf("   refit %-6s -> actual=%-5s reserve NRMSE=%.3f  ||theta||=%.2f\n",
        m, r$method, r$nrmse, sqrt(sum(r$theta^2))))
  }
}
