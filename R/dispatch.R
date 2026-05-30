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
