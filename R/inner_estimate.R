# INNER_ESTIMATE (algorithm-spec.md 4): the one step that distinguishes
# GMDH-PMM from classical GMDH.
#
#   1. LSE warm-start fit of the KG-2 partial model.
#   2. Bootstrap-stabilized cumulant diagnostics of the warm-start residuals.
#   3. Dispatch -> {LSE, PMM2, PMM3}.
#   4. If PMM, refit with EstemPMM::lm_pmm2 / lm_pmm3 (warm-started internally).
#      Forced ablation baselines can instead route to ridge-LSE / Huber / L1.
#
# The PMM estimators (Newton-Raphson, fixed/adaptive kappa) live in EstemPMM;
# this wrapper only supplies the stabilized diagnostics and the dispatch.

#' Estimate one KG-2 partial model with automatic LSE/PMM dispatch
#'
#' @param v1,v2 numeric input vectors (training partition).
#' @param y training response vector.
#' @param control a \code{\link{gmdh_pmm_control}} list.
#' @return list with \code{theta} (length-6 canonical coefficients),
#'   \code{method}, \code{diag} (cumulant diagnostics), and \code{converged}
#'   (FALSE if a PMM solver failed and the fit fell back to LSE).
#' @export
inner_estimate <- function(v1, v2, y, control = gmdh_pmm_control()) {
  d <- kg2_design(v1, v2, y)

  fit_lse <- stats::lm(.kg2_formula, data = d)
  theta <- .extract_kg2_coef(fit_lse)
  eps <- stats::residuals(fit_lse)

  diag <- bootstrap_cumulant_diag(eps, B = control$B, robust = control$robust,
                                  seed = control$boot_seed)
  method <- control$force_method %||% "auto"
  diag_weak <- NULL
  if (method == "auto") {
    method <- dispatch_method(diag, alpha = control$alpha,
                              skew_min = control$skew_min,
                              skew_strong = control$skew_strong,
                              g2_threshold = control$g2_threshold,
                              kurt_threshold = control$kurt_threshold)
  } else if (method == "auto-weak") {
    diag_weak <- weak_cumulant_diag(eps, B = control$B,
                                    sigma_mult = control$weak_sigma_mult %||% 2.5,
                                    robust = control$robust, seed = control$boot_seed)
    method <- dispatch_method_weak(diag_weak, diag, alpha = control$alpha,
                                   skew_min = control$skew_min,
                                   skew_strong = control$skew_strong,
                                   g2_threshold = control$g2_threshold,
                                   kurt_threshold = control$kurt_threshold,
                                   kurt_heavy_skew = control$kurt_heavy_skew %||% 8.0,
                                   kurt_heavy_sym = control$kurt_heavy_sym %||% 3.0)
  } else if (method == "auto-valgate") {
    # Validation-gated dispatch: pick the inner estimator by held-out error, then
    # (below) refit it on the full node training set. No cumulant thresholds -- the
    # GMDH selection principle applied to the inner estimator (a discrete super-learner).
    # Inner CV is k-fold over CONTIGUOUS blocks: this respects autocorrelation on
    # time-series nodes (a single random sub-split leaks and mis-selects there, the
    # 2026-05-31 finding) and is harmless for cross-sectional nodes (row order
    # exchangeable); averaging over k folds also cuts the selection variance of the
    # old single 70/30 split.
    cands <- control$valgate_candidates %||% c("LSE", "Huber", "L1", "WPMM2", "WPMM3")
    n <- length(y)
    method <- "LSE"
    kf <- as.integer(control$valgate_folds %||% 4L)
    if (n >= 40L && kf >= 2L) {
      .trm <- function(e) { e <- e[is.finite(e)]; if (!length(e)) return(Inf)
        q <- stats::quantile(abs(e), 0.9, names = FALSE); sqrt(mean(e[abs(e) <= q]^2)) }
      sub_ctrl <- control; sub_ctrl$B <- 0L
      fb <- n %/% kf
      # "random" is the inner-fold ablation: the same fold sizes over a fixed
      # pseudo-random permutation, so temporally adjacent rows fall on both sides.
      ord <- if (identical(control$valgate_inner, "random")) {
        order((sin(seq_len(n) * 12.9898) * 43758.5453) %% 1)
      } else seq_len(n)
      err <- matrix(NA_real_, kf, length(cands), dimnames = list(NULL, cands))
      for (k in seq_len(kf)) {
        a <- (k - 1L) * fb + 1L; z <- if (k == kf) n else k * fb
        sv <- ord[a:z]; si <- setdiff(seq_len(n), sv)            # held-out block
        if (length(si) < 12L || length(sv) < 4L) next
        for (j in seq_along(cands)) {
          sub_ctrl$force_method <- cands[j]
          th <- tryCatch(inner_estimate(v1[si], v2[si], y[si], sub_ctrl)$theta,
                         error = function(e) rep(NA_real_, 6))
          if (!any(!is.finite(th))) err[k, j] <- .trm(y[sv] - kg2_predict(th, v1[sv], v2[sv]))
        }
      }
      ms <- apply(err, 2, function(c) if (all(is.na(c))) Inf else mean(c, na.rm = TRUE))
      if (any(is.finite(ms))) method <- cands[which.min(ms)]
    }
  }
  converged <- TRUE

  if (method == "ridge-LSE") {
    theta <- .fit_kg2_ridge(d, lambda = control$ridge_lambda %||% 1e-8)

  } else if (method == "Huber") {
    theta <- .fit_kg2_irls(d, loss = "Huber", huber_k = control$huber_k %||% 1.345,
                           max_iter = control$max_iter, tol = control$tol,
                           lambda = control$ridge_lambda %||% 1e-8)

  } else if (method == "L1") {
    theta <- .fit_kg2_irls(d, loss = "L1", max_iter = control$max_iter,
                           tol = control$tol, lambda = control$ridge_lambda %||% 1e-8)

  } else if (method == "WPMM2") {
    theta <- tryCatch(
      fit_kg2_weak_pmm2(d, sigma_mult = control$weak_sigma_mult %||% 2.5,
                        max_iter = control$max_iter, tol = control$tol),
      error = function(e) NULL)
    if (is.null(theta)) { theta <- .extract_kg2_coef(fit_lse); method <- "LSE"; converged <- FALSE }

  } else if (method == "WPMM3") {
    theta <- tryCatch(
      fit_kg2_weak_pmm3(d, sigma_mult = control$weak_sigma_mult %||% 2.5,
                        max_iter = control$max_iter, tol = control$tol),
      error = function(e) NULL)
    if (is.null(theta)) { theta <- .extract_kg2_coef(fit_lse); method <- "LSE"; converged <- FALSE }

  } else if (method == "PATP3") {
    theta <- tryCatch(
      fit_kg2_patp3(d, alpha = control$patp_alpha, max_iter = control$max_iter, tol = control$tol),
      error = function(e) NULL)
    if (is.null(theta)) { theta <- .extract_kg2_coef(fit_lse); method <- "LSE"; converged <- FALSE }

  } else if (method == "PMM2") {
    fit <- tryCatch(
      EstemPMM::lm_pmm2(.kg2_formula, data = d,
                        max_iter = control$max_iter, tol = control$tol,
                        verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit)) theta <- .extract_kg2_coef(fit)
    else { method <- "LSE"; converged <- FALSE }    # robust fallback in the cascade

  } else if (method == "PMM3") {
    fit <- tryCatch(
      EstemPMM::lm_pmm3(.kg2_formula, data = d,
                        max_iter = control$max_iter, tol = control$tol,
                        adaptive = identical(control$pmm_mode, "adaptive"),
                        verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit)) theta <- .extract_kg2_coef(fit)
    else { method <- "LSE"; converged <- FALSE }
  }

  list(theta = theta, method = method, diag = diag, diag_weak = diag_weak,
       converged = converged)
}
