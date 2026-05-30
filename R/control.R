# Tuning parameters for the GMDH-PMM tournament (defaults from algorithm-spec 9).

#' GMDH-PMM control parameters
#'
#' @param split_ratio fraction of the sample used for training (default 0.6).
#' @param F top-F partial models carried to the next layer (default 10).
#' @param L_max maximum tournament depth (default 10).
#' @param epsilon relative-improvement stopping tolerance (default 1e-3): the
#'   tournament stops when the best external criterion improves by less than
#'   this fraction of the previous best.
#' @param alpha dispatch significance level (default 0.05).
#' @param B bootstrap resamples for cumulant stabilization (default 500).
#' @param max_iter Newton-Raphson iterations for PMM solvers (default 50).
#' @param tol Newton-Raphson convergence tolerance (default 1e-6).
#' @param criterion external criterion: "MSE" (default), "PMM-loss",
#'   "reserve-aware" (actuarial C5 experiments), or "spike-aware" (industrial
#'   emissions C5 experiments).
#' @param pmm_mode PMM3 kappa handling: "fixed" (default) or "adaptive".
#' @param force_method estimation override: "auto" (default), "LSE",
#'   "ridge-LSE", "Huber", "L1", "PMM2", or "PMM3". Used for ablations
#'   and baselines.
#' @param ridge_lambda ridge penalty for "ridge-LSE" (default 1e-8, matching
#'   EstemPMM's PMM2 regularization scale). The intercept is not penalized.
#' @param huber_k Huber tuning constant for the "Huber" baseline (default 1.345).
#' @param reserve_weight aggregate reserve penalty for \code{criterion =
#'   "reserve-aware"} (default 1).
#' @param tail_weight tail reserve penalty for \code{criterion =
#'   "reserve-aware"} (default 2).
#' @param overreserve_weight multiplier for over-reserve penalties relative to
#'   under-reserve penalties (default 0.25).
#' @param tail_prob quantile defining the tail subset for reserve-aware
#'   selection (default 0.9).
#' @param train_index optional integer row indices for a fixed internal training
#'   partition. If supplied, \code{val_index} must also be supplied. Used for
#'   time-aware validation experiments.
#' @param val_index optional integer row indices for a fixed internal validation
#'   partition. If supplied, \code{train_index} must also be supplied.
#' @param robust use robust (median/MAD) bootstrap summary (default TRUE).
#' @param skew_min minimum |gamma3| to consider PMM2 (default 0.3).
#' @param skew_strong |gamma3| above which PMM2 is taken outright (default 1.0).
#' @param g2_threshold g2 cutoff for PMM2 (default 0.95).
#' @param kurt_threshold gamma4 cutoff for PMM3 (default -0.7).
#' @param seed RNG seed for the train/validation split (default NULL).
#' @param boot_seed RNG seed passed to the bootstrap diagnostics (default NULL).
#' @return a named list of validated control parameters.
#' @export
gmdh_pmm_control <- function(split_ratio = 0.6, F = 10L, L_max = 10L,
                             epsilon = 1e-3, alpha = 0.05, B = 500L,
                             max_iter = 50L, tol = 1e-6,
                             criterion = c("MSE", "PMM-loss", "reserve-aware", "spike-aware"),
                             pmm_mode = c("fixed", "adaptive"),
                             force_method = c("auto", "LSE", "ridge-LSE", "Huber", "L1", "PMM2", "PMM3", "WPMM2"),
                             weak_sigma_mult = 2.5,
                             ridge_lambda = 1e-8, huber_k = 1.345,
                             reserve_weight = 1, tail_weight = 2,
                             overreserve_weight = 0.25, tail_prob = 0.9,
                             train_index = NULL, val_index = NULL,
                             robust = TRUE, skew_min = 0.3, skew_strong = 1.0,
                             g2_threshold = 0.95, kurt_threshold = -0.7,
                             seed = NULL, boot_seed = NULL) {
  criterion <- match.arg(criterion)
  pmm_mode <- match.arg(pmm_mode)
  force_method <- match.arg(force_method)
  stopifnot(split_ratio > 0, split_ratio < 1, F >= 1, L_max >= 1, B >= 0)
  stopifnot(is.finite(ridge_lambda), ridge_lambda >= 0)
  stopifnot(is.finite(huber_k), huber_k > 0)
  stopifnot(is.finite(reserve_weight), reserve_weight >= 0)
  stopifnot(is.finite(tail_weight), tail_weight >= 0)
  stopifnot(is.finite(overreserve_weight), overreserve_weight >= 0)
  stopifnot(is.finite(tail_prob), tail_prob > 0, tail_prob < 1)
  if (xor(is.null(train_index), is.null(val_index))) {
    stop("train_index and val_index must be supplied together.")
  }
  if (!is.null(train_index)) {
    stopifnot(length(train_index) >= 6, length(val_index) >= 1)
    stopifnot(all(is.finite(train_index)), all(is.finite(val_index)))
  }
  list(
    split_ratio = split_ratio, F = as.integer(F), L_max = as.integer(L_max),
    epsilon = epsilon, alpha = alpha, B = as.integer(B),
    max_iter = as.integer(max_iter), tol = tol,
    criterion = criterion, pmm_mode = pmm_mode,
    force_method = force_method, ridge_lambda = ridge_lambda,
    huber_k = huber_k, reserve_weight = reserve_weight,
    tail_weight = tail_weight, overreserve_weight = overreserve_weight,
    tail_prob = tail_prob, train_index = train_index, val_index = val_index,
    robust = robust, weak_sigma_mult = weak_sigma_mult,
    skew_min = skew_min, skew_strong = skew_strong,
    g2_threshold = g2_threshold, kurt_threshold = kurt_threshold,
    seed = seed, boot_seed = boot_seed
  )
}
