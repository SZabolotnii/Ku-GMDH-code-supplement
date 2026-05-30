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
  if (method == "auto") {
    method <- dispatch_method(diag, alpha = control$alpha,
                              skew_min = control$skew_min,
                              skew_strong = control$skew_strong,
                              g2_threshold = control$g2_threshold,
                              kurt_threshold = control$kurt_threshold)
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

  list(theta = theta, method = method, diag = diag, converged = converged)
}
