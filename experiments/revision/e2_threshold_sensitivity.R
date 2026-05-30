#!/usr/bin/env Rscript
# EXPERIMENT E2 -- Dispatch threshold sensitivity (reviewer Q1).
#
# For each of the 7 residual regimes from run_synthetic.R, generate n=2000
# residuals over 30 seeds, compute bootstrap_cumulant_diag (B=200), and record
# the dispatched method via dispatch_method() under the DEFAULT thresholds.
# Then sweep, one parameter at a time:
#   g2_threshold  in {0.90, 0.95, 0.99}
#   kurt_threshold in {-0.5, -0.7, -0.9}
#   skew_min      in {0.2, 0.3, 0.4}
#   alpha         in {0.01, 0.05, 0.10}
# For each swept value we report, per regime, the MODAL dispatched method over
# the 30 seeds, and FLAG any regime whose modal estimator CHANGES from default.
#
# Headline claim under test: the dispatched estimator is STABLE across these
# ranges (only borderline regimes may flip).
#
# Run from the package dir:
#   Rscript experiments/revision/e2_threshold_sensitivity.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))

N_RESID <- 2000L
N_SEED  <- 30L
B_BOOT  <- 200L

# --- standardized noise generators (mean 0, variance 1), identical to
#     run_synthetic.R ---------------------------------------------------------
std <- function(x) (x - mean(x)) / stats::sd(x)
gens <- list(
  G        = function(n) rnorm(n),                                   # Gaussian
  exp      = function(n) rexp(n) - 1,                                # strong + skew
  lognorm  = function(n) std(rlnorm(n, 0, 0.5)),                     # A3-like skew
  uniform  = function(n) runif(n, -sqrt(3), sqrt(3)),               # B3 platykurtic
  bimodal  = function(n) std(sample(c(-1, 1), n, TRUE) * 1.3 + rnorm(n, 0, 0.45)),
  laplace  = function(n) (rexp(n) - rexp(n)) / sqrt(2),             # B1 leptokurtic (decoy)
  contam   = function(n) std(ifelse(runif(n) < 0.9, rnorm(n), rnorm(n, 0, 5)))  # C1 (decoy)
)

# --- precompute one bootstrap diagnostic per (regime, seed); the sweep only
#     changes dispatch_method() thresholds, NOT the diagnostics, so we can
#     reuse the same diag objects across all swept values --------------------
cat(sprintf("=== E2 dispatch threshold sensitivity (n=%d, seeds=%d, B=%d) ===\n\n",
            N_RESID, N_SEED, B_BOOT))
t0 <- Sys.time()

diags <- list()  # diags[[regime]] = list of N_SEED diag objects
for (rg in names(gens)) {
  gen <- gens[[rg]]
  dl <- vector("list", N_SEED)
  for (s in seq_len(N_SEED)) {
    set.seed(7000L + s)          # seed for residual generation
    e <- gen(N_RESID)
    # boot seed kept distinct & deterministic per (regime, seed)
    dl[[s]] <- bootstrap_cumulant_diag(e, B = B_BOOT, robust = TRUE,
                                       seed = 9000L + s)
  }
  diags[[rg]] <- dl
}
cat(sprintf("Bootstrap diagnostics computed (%.1fs).\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# helper: modal method over the 30 seeds for a given threshold config --------
modal_method <- function(dl, alpha, skew_min, skew_strong,
                         g2_threshold, kurt_threshold) {
  m <- vapply(dl, function(d)
    dispatch_method(d, alpha = alpha, skew_min = skew_min,
                    skew_strong = skew_strong, g2_threshold = g2_threshold,
                    kurt_threshold = kurt_threshold),
    character(1))
  tab <- sort(table(factor(m, levels = c("LSE", "PMM2", "PMM3"))),
              decreasing = TRUE)
  list(mode = names(tab)[1], counts = m)
}

# default thresholds ---------------------------------------------------------
DEF <- list(alpha = 0.05, skew_min = 0.3, skew_strong = 1.0,
            g2_threshold = 0.95, kurt_threshold = -0.7)

# default dispatch table + mean diagnostics (context) ------------------------
def_mode <- list()
def_rows <- list()
for (rg in names(gens)) {
  dl <- diags[[rg]]
  mm <- modal_method(dl, DEF$alpha, DEF$skew_min, DEF$skew_strong,
                     DEF$g2_threshold, DEF$kurt_threshold)
  def_mode[[rg]] <- mm$mode
  g3v <- vapply(dl, `[[`, numeric(1), "gamma3")
  g4v <- vapply(dl, `[[`, numeric(1), "gamma4")
  g2v <- vapply(dl, `[[`, numeric(1), "g2")
  tab <- table(factor(mm$counts, levels = c("LSE", "PMM2", "PMM3")))
  def_rows[[rg]] <- data.frame(
    regime      = rg,
    mean_gamma3 = round(mean(g3v), 3),
    mean_gamma4 = round(mean(g4v), 3),
    mean_g2     = round(mean(g2v), 3),
    n_LSE       = as.integer(tab["LSE"]),
    n_PMM2      = as.integer(tab["PMM2"]),
    n_PMM3      = as.integer(tab["PMM3"]),
    modal       = mm$mode,
    stringsAsFactors = FALSE
  )
}
def_table <- do.call(rbind, def_rows)
cat("DEFAULT dispatch table (thresholds: alpha=0.05, skew_min=0.3,",
    "skew_strong=1.0, g2_threshold=0.95, kurt_threshold=-0.7):\n")
print(def_table, row.names = FALSE)
cat("\n")

# --- one-at-a-time sweep ----------------------------------------------------
sweeps <- list(
  g2_threshold   = c(0.90, 0.95, 0.99),
  kurt_threshold = c(-0.5, -0.7, -0.9),
  skew_min       = c(0.2, 0.3, 0.4),
  alpha          = c(0.01, 0.05, 0.10)
)

rows <- list()
flag_lines <- character(0)
for (param in names(sweeps)) {
  for (val in sweeps[[param]]) {
    cfg <- DEF
    cfg[[param]] <- val
    for (rg in names(gens)) {
      mm <- modal_method(diags[[rg]], cfg$alpha, cfg$skew_min, cfg$skew_strong,
                         cfg$g2_threshold, cfg$kurt_threshold)
      changed <- mm$mode != def_mode[[rg]]
      tab <- table(factor(mm$counts, levels = c("LSE", "PMM2", "PMM3")))
      rows[[length(rows) + 1L]] <- data.frame(
        param         = param,
        value         = val,
        regime        = rg,
        modal_default = def_mode[[rg]],
        modal_swept   = mm$mode,
        changed       = changed,
        n_LSE         = as.integer(tab["LSE"]),
        n_PMM2        = as.integer(tab["PMM2"]),
        n_PMM3        = as.integer(tab["PMM3"]),
        stringsAsFactors = FALSE
      )
      if (changed) {
        flag_lines <- c(flag_lines, sprintf(
          "  FLAG: %s=%s  regime=%s  %s -> %s",
          param, format(val), rg, def_mode[[rg]], mm$mode))
      }
    }
  }
}
sweep_df <- do.call(rbind, rows)

# --- per-parameter summary --------------------------------------------------
for (param in names(sweeps)) {
  cat(sprintf("--- Sweep %s in {%s} ---\n", param,
              paste(format(sweeps[[param]]), collapse = ", ")))
  sub <- sweep_df[sweep_df$param == param, ]
  # wide table: regime x value -> modal method, '*' marks a change
  vals <- sweeps[[param]]
  hdr <- sprintf("%-9s", "regime")
  for (v in vals) hdr <- paste0(hdr, sprintf("%12s", format(v)))
  cat(hdr, "\n")
  for (rg in names(gens)) {
    line <- sprintf("%-9s", rg)
    for (v in vals) {
      r <- sub[sub$regime == rg & sub$value == v, ]
      mark <- if (isTRUE(r$changed)) "*" else " "
      line <- paste0(line, sprintf("%11s%s", r$modal_swept, mark))
    }
    cat(line, "\n")
  }
  chg <- sub[sub$changed, ]
  if (nrow(chg) == 0) {
    cat("  -> NO regime changes from default across this sweep.\n\n")
  } else {
    cat(sprintf("  -> %d (regime,value) change(s) from default: %s\n\n",
                nrow(chg),
                paste(unique(chg$regime), collapse = ", ")))
  }
}

# --- overall flags ----------------------------------------------------------
cat("=== FLAGGED changes (modal estimator differs from default) ===\n")
if (length(flag_lines) == 0) {
  cat("  NONE -- dispatched estimator is STABLE across all swept ranges.\n")
} else {
  cat(paste(flag_lines, collapse = "\n"), "\n")
}

# headline verdict
n_changed <- sum(sweep_df$changed)
regimes_changed <- unique(sweep_df$regime[sweep_df$changed])
cat(sprintf("\n>>> %d of %d swept (regime,value) cells change from default.\n",
            n_changed, nrow(sweep_df)))
cat(sprintf(">>> Regimes that ever flip: %s\n",
            if (length(regimes_changed) == 0) "(none)"
            else paste(regimes_changed, collapse = ", ")))

# --- persist ----------------------------------------------------------------
outdir <- "experiments/revision"
utils::write.csv(def_table, file.path(outdir, "e2_default_dispatch.csv"),
                 row.names = FALSE)
utils::write.csv(sweep_df, file.path(outdir, "e2_threshold_sensitivity.csv"),
                 row.names = FALSE)
cat(sprintf("\nSaved %s and %s | elapsed %.1fs\n",
            file.path(outdir, "e2_default_dispatch.csv"),
            file.path(outdir, "e2_threshold_sensitivity.csv"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
