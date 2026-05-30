#!/usr/bin/env Rscript
# EXPERIMENT E4 -- LEVERAGE / BOOTSTRAP ROBUSTNESS (reviewer Q3).
#
# Reviewer concern: the iid residual bootstrap behind bootstrap_cumulant_diag()
# resamples raw LSE residuals e_i, which are heteroskedastic under leverage:
# Var(e_i) = sigma^2 (1 - h_ii). High-leverage rows therefore have ARTIFICIALLY
# SMALL raw residuals, which can deflate (or otherwise distort) the estimated
# skewness gamma3 / excess kurtosis gamma4 and so flip the LSE/PMM2/PMM3
# dispatch. The classical fix is to studentize the residuals before resampling:
#   e_i^stud = e_i / sqrt(1 - h_ii),   h_ii = diag of the hat matrix H = X(X'X)^-1 X'.
# These leverage-corrected residuals are (asymptotically) homoskedastic, so an
# iid bootstrap on them is justified.
#
# DESIGN. For each of 5 representative residual regimes we build a KG-2 design
# (b1,b2,b12,b11,b22) on n points with NON-TRIVIAL leverage: a handful of rows
# are pushed far out in v1/v2, which the quadratic features b11/b22 amplify into
# very high h_ii. We fit a single LSE KG-2 model, then estimate standardized
# residual cumulants two ways, holding the bootstrap RNG seed fixed so the ONLY
# difference is raw vs studentized residuals:
#   (a) iid  : bootstrap_cumulant_diag(e)
#   (b) studz: bootstrap_cumulant_diag(e / sqrt(1 - h_ii))
# and record the dispatched method under each. Over 30 seeds we report mean
# gamma3, gamma4 and the dispatch agreement rate per regime.
#
# Usage: Rscript experiments/revision/e4_leverage_bootstrap.R [N_SEEDS] [B]

suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(".", quiet = TRUE)
)

args    <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 30L
B       <- if (length(args) >= 2) as.integer(args[2]) else 200L
N       <- 400L                                  # rows per design
SEED0   <- 91001L
THETA   <- c(0.5, 1.5, -1.0, 0.8, 0.6, -0.4)     # canonical KG-2 truth
SIGMA   <- 1.0

# --- standardized noise generators (mean 0, variance 1), mirroring run_synthetic.R
std <- function(x) (x - mean(x)) / stats::sd(x)
gens <- list(
  exp_skew      = function(n) rexp(n) - 1,                                       # A-type asymmetric -> PMM2
  lognorm_skew  = function(n) std(rlnorm(n, 0, 0.5)),                            # A3-like skew -> PMM2
  uniform_platy = function(n) runif(n, -sqrt(3), sqrt(3)),                       # B3 platykurtic -> PMM3
  bimodal_platy = function(n) std(sample(c(-1, 1), n, TRUE) * 1.3 + rnorm(n, 0, 0.45)), # platykurtic -> PMM3
  contam_decoy  = function(n) std(ifelse(runif(n) < 0.9, rnorm(n), rnorm(n, 0, 5)))     # C1 leptokurtic decoy -> LSE
)

# Build a KG-2 design with deliberate leverage: ~4% of rows are pushed far out
# along v1 or v2; the quadratic KG-2 features (v1^2, v2^2) then turn these into
# genuinely high-leverage rows in the model.matrix sense (large h_ii).
make_leverage_design <- function(n, seed) {
  set.seed(seed)
  v1 <- rnorm(n); v2 <- rnorm(n)
  k <- max(4L, round(0.02 * n))                  # high-leverage rows per axis
  i1 <- seq_len(k); i2 <- k + seq_len(k)
  v1[i1] <- v1[i1] + sample(c(-1, 1), k, TRUE) * runif(k, 6, 10)
  v2[i2] <- v2[i2] + sample(c(-1, 1), k, TRUE) * runif(k, 6, 10)
  list(v1 = v1, v2 = v2)
}

# Hat-matrix diagonals for the KG-2 lm design (intercept + 5 features, p = 6).
hat_diag <- function(v1, v2) {
  X <- stats::model.matrix(gmdhpmm:::.kg2_formula, data = kg2_design(v1, v2, y = rep(0, length(v1))))
  XtX <- crossprod(X)
  H <- X %*% solve(XtX + 1e-10 * diag(ncol(X))) %*% t(X)
  list(h = pmin(pmax(diag(H), 0), 1 - 1e-8), p = ncol(X))
}

run_cell <- function(name, gen, seed) {
  des <- make_leverage_design(N, seed)
  v1 <- des$v1; v2 <- des$v2
  set.seed(seed + 5e5)                            # noise RNG, independent of design
  y  <- kg2_predict(THETA, v1, v2) + SIGMA * gen(N)

  d   <- kg2_design(v1, v2, y)
  fit <- stats::lm(gmdhpmm:::.kg2_formula, data = d)
  e   <- stats::residuals(fit)

  hd  <- hat_diag(v1, v2)
  h   <- hd$h
  e_stud <- e / sqrt(1 - h)                        # leverage-corrected (studentized) residuals

  # Hold the bootstrap RNG seed fixed across the two variants: the ONLY thing
  # that changes is raw vs studentized residuals.
  bseed <- seed + 7777L
  d_iid <- bootstrap_cumulant_diag(e,      B = B, robust = TRUE, seed = bseed)
  d_stu <- bootstrap_cumulant_diag(e_stud, B = B, robust = TRUE, seed = bseed)

  m_iid <- dispatch_method(d_iid)
  m_stu <- dispatch_method(d_stu)

  data.frame(
    regime       = name,
    seed         = seed,
    n            = N,
    p            = hd$p,
    max_h        = max(h),
    mean_h       = mean(h),
    frac_hi_lev  = mean(h > 2 * hd$p / N),         # rule-of-thumb high-leverage fraction (h > 2p/n)
    g3_iid       = d_iid$gamma3,
    g3_stud      = d_stu$gamma3,
    g4_iid       = d_iid$gamma4,
    g4_stud      = d_stu$gamma4,
    se_g3_iid    = d_iid$se_gamma3,
    se_g3_stud   = d_stu$se_gamma3,
    method_iid   = m_iid,
    method_stud  = m_stu,
    agree        = m_iid == m_stu,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("=== E4 leverage / bootstrap robustness (seeds=%d, n=%d, B=%d) ===\n\n",
            N_SEEDS, N, B))
t0 <- Sys.time()

raw <- do.call(rbind, lapply(names(gens), function(nm)
  do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_cell(nm, gens[[nm]], SEED0 + s - 1L)))))

# --- per-regime summary ----------------------------------------------------
agg <- do.call(rbind, by(raw, raw$regime, function(g) {
  mode1 <- function(x) names(sort(table(x), decreasing = TRUE))[1]
  data.frame(
    regime        = g$regime[1],
    n_seed        = nrow(g),
    max_h         = round(max(g$max_h), 3),
    mean_frac_hi  = round(mean(g$frac_hi_lev), 4),
    g3_iid_mean   = round(mean(g$g3_iid), 3),
    g3_stud_mean  = round(mean(g$g3_stud), 3),
    g3_shift      = round(mean(g$g3_stud - g$g3_iid), 4),
    g4_iid_mean   = round(mean(g$g4_iid), 3),
    g4_stud_mean  = round(mean(g$g4_stud), 3),
    g4_shift      = round(mean(g$g4_stud - g$g4_iid), 4),
    disp_iid      = mode1(g$method_iid),
    disp_stud     = mode1(g$method_stud),
    agree_rate    = round(mean(g$agree), 3),
    stringsAsFactors = FALSE
  )
}))
agg <- agg[match(names(gens), agg$regime), ]

cat("Per-regime summary (mean over seeds; *_shift = stud - iid):\n")
print(agg[, c("regime","n_seed","max_h","mean_frac_hi",
              "g3_iid_mean","g3_stud_mean","g3_shift",
              "g4_iid_mean","g4_stud_mean","g4_shift",
              "disp_iid","disp_stud","agree_rate")],
      row.names = FALSE, digits = 4)

cat("\nDispatch decision tables (iid vs studentized), per regime:\n")
for (nm in names(gens)) {
  g <- raw[raw$regime == nm, ]
  cat(sprintf("\n  [%s]\n", nm))
  print(table(iid = g$method_iid, studentized = g$method_stud))
}

overall_agree <- mean(raw$agree)
cat(sprintf("\n>>> Overall dispatch agreement (iid vs leverage-corrected): %.1f%% (%d/%d cells)\n",
            100 * overall_agree, sum(raw$agree), nrow(raw)))

# --- persist ---------------------------------------------------------------
outdir <- "experiments/revision"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(raw, file.path(outdir, "e4_leverage_bootstrap.csv"), row.names = FALSE)
utils::write.csv(agg, file.path(outdir, "e4_leverage_bootstrap_summary.csv"), row.names = FALSE)

cat(sprintf("\nSaved %s\n      %s\n      | elapsed %.1fs\n",
            file.path(outdir, "e4_leverage_bootstrap.csv"),
            file.path(outdir, "e4_leverage_bootstrap_summary.csv"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
