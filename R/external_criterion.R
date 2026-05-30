# External selection criterion (algorithm-spec.md 8).
#
# Default is the classical MSE regularity criterion. The optional "PMM-loss"
# criterion uses the integral PMM loss matching the chosen method; it is
# experimental (8.2 flags comparability concerns across methods), so MSE stays
# the default and is always reported alongside. The "reserve-aware" criterion
# is an actuarial C5 experiment: it augments MSE with aggregate and tail reserve
# adequacy penalties on the validation set. The "spike-aware" criterion targets
# industrial-emission C5 cases: it also upweights pointwise under-prediction in
# the observed high-emission tail.

#' External regularity criterion on the validation set
#'
#' @param y_val observed validation responses.
#' @param yhat_val partial-model predictions on the validation set.
#' @param method estimation method used for the partial model
#'   ("LSE", "PMM2", "PMM3"); only matters for \code{type = "PMM-loss"}.
#' @param diag diagnostics list from \code{\link{bootstrap_cumulant_diag}}
#'   (supplies the cumulants/moments for the PMM loss).
#' @param type "MSE" (default), "PMM-loss", "reserve-aware", or
#'   "spike-aware".
#' @param reserve_weight aggregate reserve penalty weight for
#'   \code{type = "reserve-aware"}.
#' @param tail_weight tail reserve penalty weight for
#'   \code{type = "reserve-aware"}.
#' @param overreserve_weight over-reserve penalty multiplier relative to
#'   under-reserve. Values below 1 encode the actuarial preference that
#'   under-reserving is more dangerous than mild over-reserving.
#' @param tail_prob quantile defining the high-claim tail subset.
#' @return scalar criterion value (lower is better).
#' @export
external_criterion <- function(y_val, yhat_val, method = "LSE", diag = NULL,
                               type = c("MSE", "PMM-loss", "reserve-aware", "spike-aware"),
                               reserve_weight = 1, tail_weight = 2,
                               overreserve_weight = 0.25, tail_prob = 0.9) {
  type <- match.arg(type)
  e <- y_val - yhat_val
  mse <- mean(e^2)
  if (type %in% c("reserve-aware", "spike-aware")) {
    scale2 <- stats::var(y_val)
    if (!is.finite(scale2) || scale2 <= 1e-12) scale2 <- mean(y_val^2)
    if (!is.finite(scale2) || scale2 <= 1e-12) scale2 <- 1

    reserve_error <- .reserve_error(y_val, yhat_val)
    penalty <- .asymmetric_square_penalty(reserve_error, overreserve_weight)

    q <- as.numeric(stats::quantile(y_val, tail_prob, na.rm = TRUE, names = FALSE))
    tail <- is.finite(q) & y_val >= q
    tail_error <- if (any(tail, na.rm = TRUE)) .reserve_error(y_val[tail], yhat_val[tail]) else 0
    tail_penalty <- .asymmetric_square_penalty(tail_error, overreserve_weight)

    cr <- mse + scale2 * (reserve_weight * penalty + tail_weight * tail_penalty)
    if (type == "spike-aware" && any(tail, na.rm = TRUE)) {
      tail_e <- y_val[tail] - yhat_val[tail]
      under <- pmax(0, tail_e)
      over <- pmax(0, -tail_e)
      spike_penalty <- mean(under^2 + overreserve_weight * over^2)
      cr <- cr + tail_weight * spike_penalty
    }
    return(cr)
  }
  if (type == "MSE" || method == "LSE" || is.null(diag)) return(mse)

  pt <- diag$point
  if (method == "PMM2") {
    c2 <- pt$m2
    Delta <- c2 * (2 + diag$gamma4 - diag$gamma3^2)
    if (!is.finite(Delta) || abs(Delta) < 1e-12) return(mse)
    a1 <- (2 + diag$gamma4) / Delta
    a2 <- -diag$gamma3 * sqrt(c2) / Delta
    return(mean(a1 * e^2 / 2 + a2 * e^3 / 3))
  }
  if (method == "PMM3") {
    denom <- pt$m4 - 3 * pt$m2^2
    if (!is.finite(denom) || abs(denom) < 1e-12) return(mse)
    kappa <- (pt$m6 - 3 * pt$m4 * pt$m2) / denom
    return(mean(kappa * e^2 / 2 - e^4 / 4))
  }
  mse
}

.reserve_error <- function(y, yhat) {
  den <- sum(y)
  if (!is.finite(den) || abs(den) < 1e-12) return(0)
  out <- sum(yhat) / den - 1
  if (is.finite(out)) out else 0
}

.asymmetric_square_penalty <- function(reserve_error, overreserve_weight) {
  if (!is.finite(reserve_error)) return(0)
  under <- max(0, -reserve_error)
  over <- max(0, reserve_error)
  under^2 + overreserve_weight * over^2
}
