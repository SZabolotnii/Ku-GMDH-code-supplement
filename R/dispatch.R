# Automatic LSE / PMM2 / PMM3 dispatch (algorithm-spec.md 4.1).
#
# Difference from EstemPMM::pmm_dispatch(): this version gates on statistical
# significance of the bootstrap-stabilized cumulants (|gamma3| > z * SE), which
# suppresses false positives from sampling noise on small GMDH partitions.

#' Select an estimation method from cumulant diagnostics
#'
#' Implements the decision flow of algorithm-spec.md 4.1:
#' \itemize{
#'   \item significant skew (\eqn{|\gamma_3| > z\,SE} and \eqn{|\gamma_3| > 0.3})
#'     with \eqn{|\gamma_3| > 1.0} or \eqn{g_2 < 0.95} -> \code{"PMM2"};
#'   \item else symmetric and significantly platykurtic
#'     (\eqn{\gamma_4 < -0.7} and \eqn{|\gamma_4| > z\,SE}) -> \code{"PMM3"};
#'   \item otherwise -> \code{"LSE"}.
#' }
#'
#' @param diag a diagnostics list from \code{\link{bootstrap_cumulant_diag}}.
#' @param alpha dispatch significance level (default 0.05).
#' @param skew_min minimum |gamma3| to consider PMM2 (default 0.3).
#' @param skew_strong |gamma3| above which PMM2 is taken outright (default 1.0).
#' @param g2_threshold g2 below which PMM2 is worthwhile (default 0.95).
#' @param kurt_threshold gamma4 below which PMM3 applies (default -0.7).
#' @return one of \code{"LSE"}, \code{"PMM2"}, \code{"PMM3"}.
#' @export
dispatch_method <- function(diag, alpha = 0.05,
                            skew_min = 0.3, skew_strong = 1.0,
                            g2_threshold = 0.95, kurt_threshold = -0.7) {
  z <- stats::qnorm(1 - alpha / 2)
  g3 <- diag$gamma3; g4 <- diag$gamma4
  se3 <- diag$se_gamma3; se4 <- diag$se_gamma4

  skew_significant <- isTRUE(abs(g3) > z * se3) && isTRUE(abs(g3) > skew_min)
  if (skew_significant) {
    if (isTRUE(abs(g3) > skew_strong) || isTRUE(diag$g2 < g2_threshold)) return("PMM2")
    return("LSE")
  }

  kurt_significant <- isTRUE(g4 < kurt_threshold) && isTRUE(abs(g4) > z * se4)
  if (kurt_significant) return("PMM3")

  "LSE"
}

#' Select a WEAK estimation method from weak + raw cumulant diagnostics
#'
#' The weak dispatch (Paper 4). Weak-polynomial corrections (windowed PMM) only
#' help when there is a genuine heavy tail for the window to tame; a regime map
#' showed the weak skew correction (WPMM2) is harmful under light-tailed bulk skew
#' and inert under symmetry, while the weak kurtosis correction (WPMM3) gives
#' bandwidth-insurance under symmetric heavy tails. The decision flow is therefore:
#' \itemize{
#'   \item heavy tail (\code{tail_ratio > tail_ratio_heavy} or
#'     \eqn{\gamma_4^{raw} > kurt\_heavy}) with significant weak skew
#'     (\eqn{|\gamma_3^w| > z\,SE} and \eqn{|\gamma_3^w| > skew\_min}, and
#'     \eqn{|\gamma_3^w| > skew\_strong} or \eqn{g_2^w < g2\_threshold})
#'     -> \code{"WPMM2"};
#'   \item heavy tail, symmetric -> \code{"WPMM3"} (bandwidth insurance);
#'   \item light tail -> defer to the classical \code{\link{dispatch_method}}
#'     (LSE / classical PMM2 on clean skew / raw PMM3 on platykurtic).
#' }
#'
#' @param diag_weak weak diagnostics from \code{\link{weak_cumulant_diag}}.
#' @param diag_raw raw diagnostics from \code{\link{bootstrap_cumulant_diag}}
#'   (used for the light-tail / clean-skew classical fallback).
#' @param alpha,skew_min,skew_strong,g2_threshold,kurt_threshold as in
#'   \code{\link{dispatch_method}} (applied to the weak skew test / classical fallback).
#' @param kurt_heavy_skew raw excess kurtosis above which a SKEWED residual is
#'   treated as asymmetric heavy contamination -> WPMM2 (default 8.0). Below it,
#'   skew is treated as intrinsic / clean -> classical PMM2 (the regime where the
#'   weak skew tilt is harmful, e.g. exponential / lognormal errors).
#' @param kurt_heavy_sym raw excess kurtosis above which a SYMMETRIC residual is
#'   treated as heavy-tailed -> WPMM3 (default 3.0). Below it -> classical
#'   (LSE, or raw PMM3 on platykurtic errors).
#' @return one of \code{"WPMM2"}, \code{"WPMM3"}, \code{"LSE"}, \code{"PMM2"}, \code{"PMM3"}.
#' @export
dispatch_method_weak <- function(diag_weak, diag_raw, alpha = 0.05,
                                 skew_min = 0.3, skew_strong = 1.0,
                                 g2_threshold = 0.95, kurt_threshold = -0.7,
                                 kurt_heavy_skew = 8.0, kurt_heavy_sym = 3.0) {
  classical <- function() dispatch_method(diag_raw, alpha = alpha, skew_min = skew_min,
                                          skew_strong = skew_strong, g2_threshold = g2_threshold,
                                          kurt_threshold = kurt_threshold)
  g4_raw <- diag_weak$gamma4_raw
  if (!isTRUE(is.finite(g4_raw))) return(classical())

  z <- stats::qnorm(1 - alpha / 2)
  g3w <- diag_weak$gamma3; se3w <- diag_weak$se_gamma3
  skew_significant <- isTRUE(abs(g3w) > z * se3w) && isTRUE(abs(g3w) > skew_min)

  if (skew_significant) {
    # Skewed: the weak skew tilt (WPMM2) helps only when the asymmetry rides on a
    # genuine heavy tail / contamination; under intrinsic light-tail skew it is
    # harmful, so defer to the classical (raw) PMM2 dispatch instead.
    if (isTRUE(g4_raw > kurt_heavy_skew)) return("WPMM2")
    return(classical())
  }
  # Symmetric: WPMM3 (bandwidth insurance) under heavy tails; classical otherwise
  # (LSE on light tails, raw PMM3 on platykurtic).
  if (isTRUE(g4_raw > kurt_heavy_sym)) return("WPMM3")
  classical()
}
