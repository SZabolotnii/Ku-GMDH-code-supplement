#!/usr/bin/env Rscript
# Volve in-distribution EFFICIENCY experiment (Paper 1 revision, follow-up).
#
# The ReserveTW4 audits test EXTRAPOLATION past a formation break, where
# robustness/structure dominate and PMM is only intermediate. PMM's home turf is
# EFFICIENCY: lower-variance estimation when residuals carry stable, transferable
# non-Gaussian structure. EDA shows Volve LAS file 3 (~3618-4089 m) has strongly,
# stably skewed log-ROP residuals (gamma3 ~ -1.4, g2 ~ 0.68) that transfer from
# train to a later-depth, same-regime test window -- the setting where PMM2
# should beat LSE and ideally also LAD/Huber on a standard temporal holdout.
#
# Usage: Rscript experiments/run_volve_efficiency.R [N_SEEDS] [FILE] [SUBSAMPLE]
#   FILE = LAS index 1..4 (default 3).

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS   <- if (length(args) >= 1) as.integer(args[1]) else 30L
FILEIDX   <- if (length(args) >= 2) as.integer(args[2]) else 3L
SUBSAMPLE <- if (length(args) >= 3) as.integer(args[3]) else 0L
ROWMAX    <- if (length(args) >= 4) as.integer(args[4]) else 0L  # restrict to first ROWMAX rows
ROWMIN    <- if (length(args) >= 5) as.integer(args[5]) else 1L  # start row (heavy-tail window)
SEED0 <- 75001L
L_MAX <- 3L; F_KEEP <- 6L; LEVEL <- 0.95; B <- 200L
METHODS <- c("PMM","WPMM2","LSE","ridge-LSE","Huber","L1")

d <- fread("data/volve_onbottom.csv")
fn <- sprintf("WL_RAW_BHPR-GR-MECH_TIME_MWD_%d.LAS", FILEIDX)
seg <- d[source_file == fn][order(dept)]
hi <- if (ROWMAX > 0 && ROWMAX <= nrow(seg)) ROWMAX else nrow(seg)
seg <- seg[ROWMIN:hi]
X_all <- as.matrix(seg[, c("swob","rpm","tqa","dept")]); storage.mode(X_all) <- "double"
y_all <- log(pmax(seg$rop5, 1e-6))
N <- nrow(seg)

std <- function(X, ref){mu<-colMeans(X[ref,,drop=FALSE]);sg<-apply(X[ref,,drop=FALSE],2,sd)
  sg[!is.finite(sg)|sg<=1e-12]<-1; sweep(sweep(X,2,mu,"-"),2,sg,"/")}
nrmse <- function(y,p){s<-sd(y);if(!is.finite(s)||s<=0)s<-1;sqrt(mean((p-y)^2))/s}
ctrl_for <- function(m,n_tr,n_val,seed){
  base<-list(L_max=L_MAX,F=F_KEEP,epsilon=-Inf,train_index=seq_len(n_tr),
    val_index=n_tr+seq_len(n_val),criterion="MSE",max_iter=60L)
  switch(m,
    PMM=do.call(gmdh_pmm_control,c(base,list(B=B,force_method="auto",boot_seed=seed))),
    WPMM2=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="WPMM2",weak_sigma_mult=2.5))),
    LSE=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="LSE"))),
    `ridge-LSE`=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="ridge-LSE",ridge_lambda=1e-2))),
    Huber=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="Huber"))),
    L1=do.call(gmdh_pmm_control,c(base,list(B=0,force_method="L1"))))
}
mshare <- function(fit,m){tabs<-lapply(fit$layers,function(z)z$methods)
  num<-sum(vapply(tabs,function(t){v<-unname(t[m]);if(!length(v)||is.na(v))0 else v},numeric(1)))
  den<-sum(vapply(tabs,sum,numeric(1)));if(den<=0)NA_real_ else num/den}

# In-distribution temporal split within the stationary segment:
# 55% train, 20% selection-val, 10% interval-calibration, 15% test (later depth).
run_one <- function(seed){
  set.seed(seed)
  fit_end <- floor(0.75*N); cal_end <- floor(0.85*N)
  fit_block <- seq_len(fit_end)
  if (SUBSAMPLE>0 && SUBSAMPLE<length(fit_block)) fit_block <- sort(sample(fit_block,SUBSAMPLE))
  T <- length(fit_block); n_tr<-floor(0.73*T); n_val<-T-n_tr
  cal <- (fit_end+1L):cal_end; test <- (cal_end+1L):N
  Xs <- std(X_all, fit_block)
  rows<-list()
  for(m in METHODS){
    fit<-gmdh_pmm(Xs[fit_block,,drop=FALSE], y_all[fit_block], ctrl_for(m,n_tr,n_val,seed))
    pcal<-predict(fit,Xs[cal,,drop=FALSE]); pte<-predict(fit,Xs[test,,drop=FALSE])
    q<-as.numeric(quantile(abs(y_all[cal]-pcal),LEVEL,names=FALSE))
    rows[[length(rows)+1L]]<-data.frame(seed=seed,method=m,
      nrmse_test=nrmse(y_all[test],pte), coverage=mean(y_all[test]>=pte-q & y_all[test]<=pte+q),
      width=2*q, pmm2_share=mshare(fit,"PMM2"), pmm3_share=mshare(fit,"PMM3"),
      stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}

cat(sprintf("=== Volve EFFICIENCY (in-distribution) | LAS file %d, N=%d, seeds=%d ===\n", FILEIDX, N, N_SEEDS))
t0<-Sys.time()
all<-do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_one(SEED0+s-1L)))
outdir<-"experiments/results"; if(!dir.exists(outdir))dir.create(outdir,recursive=TRUE)
utils::write.csv(all, file.path(outdir,"volve_efficiency_raw.csv"), row.names=FALSE)
agg<-do.call(rbind, by(all, all$method, function(g)
  data.frame(method=g$method[1],n_seed=nrow(g),
    nrmse_med=median(g$nrmse_test),nrmse_q1=quantile(g$nrmse_test,.25,names=FALSE),
    nrmse_q3=quantile(g$nrmse_test,.75,names=FALSE),
    cov_med=median(g$coverage),width_med=median(g$width),
    pmm2_med=median(g$pmm2_share,na.rm=TRUE),pmm3_med=median(g$pmm3_share,na.rm=TRUE),
    stringsAsFactors=FALSE)))
agg<-agg[order(agg$nrmse_med),]
utils::write.csv(agg, file.path(outdir,"volve_efficiency_summary.csv"), row.names=FALSE)
cat("\n=== In-distribution test (median across seeds) ===\n")
print(agg, row.names=FALSE, digits=4)
cat(sprintf("\nSaved volve_efficiency_{raw,summary}.csv | %.0fs\n", as.numeric(difftime(Sys.time(),t0,"secs"))))
