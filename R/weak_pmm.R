# Weak-moment PMM2 inner estimator (proof-of-concept bridge to the
# Ku_Weak_Moment program). The classical PMM2 score uses RAW residual cumulants
# (gamma3, gamma4), which a single heavy-tail outlier destabilizes -- exactly
# the regime where, on real drilling data, PMM2 loses to Huber/LAD. Weak moments
# replace raw averaging E[g(r)] by a value-domain windowed functional
#   E_w[g(r)] = sum_i w_sigma(r_i) g(r_i) / sum_i w_sigma(r_i),
# with a Gaussian window centred on a ROBUST residual location. This keeps the
# skewness-exploiting PMM2 correction while taming heavy tails.
#
# Construction (consistent with cumulants.R / external_criterion.R):
#   window     w_i = exp(-0.5 ((r_i - mu_w) / (sigma_mult * s))^2),  s = MAD/0.6745
#   weak moms  c2^w, c3^w, c4^w (central, around the weak mean)
#   PMM2 coefs a1 = (2 + g4^w)/Delta,  a2 = -g3^w sqrt(c2^w)/Delta,
#              Delta = c2^w (2 + g4^w - (g3^w)^2)
#   score      U(theta) = Z' [ w (a1 r_c + a2 r_c^2) ] = 0   (r_c centred)
# solved by Gauss-Newton warm-started from LSE; weights are re-evaluated each
# step (a redescending, window-parametrized estimator). Degree-1 (a2 = 0) would
# reduce to a Welsch M-estimator; the a2 term is the skewness-exploiting part.

#' Fit a KG-2 partial model by weak-moment PMM2
#'
#' @param d data.frame with KG-2 features (b1,b2,b12,b11,b22) and response y.
#' @param sigma_mult window scale as a multiple of the robust residual scale
#'   (default 2.5; wide -> classical PMM2, narrow -> more robust).
#' @param max_iter,tol Gauss-Newton controls.
#' @return length-6 coefficient vector in canonical order.
#' @export
fit_kg2_weak_pmm2 <- function(d, sigma_mult = 2.5, max_iter = 50L, tol = 1e-7) {
  Z <- stats::model.matrix(.kg2_formula, data = d)
  y <- d$y
  theta <- .fit_kg2_ridge(d, lambda = 1e-8)        # LSE warm start
  # IRLS on the factored weak-PMM2 score: sum_i w_i (1 + a2 e_ic) e_i z_i = 0,
  # i.e. weighted least squares with weight W_i = w_i (1 + a2 e_ic). The window
  # w_i gives redescending robustness; the (1 + a2 e_ic) factor is the
  # skewness-exploiting PMM2 correction from weak cumulants.
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - Z %*% theta)
    mu <- stats::median(r)
    s <- stats::mad(r, constant = 1.4826)
    if (!is.finite(s) || s <= 1e-12) s <- stats::sd(r)
    if (!is.finite(s) || s <= 1e-12) s <- 1
    sw <- sigma_mult * s
    w <- exp(-0.5 * ((r - mu) / sw)^2)             # value-domain window, robust centre
    sw_sum <- sum(w)
    if (!is.finite(sw_sum) || sw_sum <= 0) break
    rc <- r - sum(w * r) / sw_sum                  # centre on weak mean
    c2 <- sum(w * rc^2) / sw_sum
    c3 <- sum(w * rc^3) / sw_sum
    c4 <- sum(w * rc^4) / sw_sum - 3 * c2^2        # weak 4th cumulant
    denom <- 2 * c2^2 + c4
    a2 <- if (!is.finite(c2) || c2 <= 1e-12 || !is.finite(denom) || abs(denom) < 1e-12)
            0 else -c3 / denom
    a2 <- max(min(a2, 0.5 / sqrt(c2 + 1e-12)), -0.5 / sqrt(c2 + 1e-12))  # keep a2*rc bounded
    W <- w * (1 + a2 * rc)
    W[!is.finite(W) | W < 0] <- 0                  # redescending: never negative weight
    theta_new <- .solve_weighted_ridge(Z, y, W, lambda = 1e-8)
    if (any(!is.finite(theta_new))) break
    if (sqrt(sum((theta_new - theta)^2)) / max(1, sqrt(sum(theta^2))) < tol) {
      theta <- theta_new; break
    }
    theta <- theta_new
  }
  out <- stats::setNames(as.numeric(theta), .kg2_coef_order)
  out[!is.finite(out)] <- 0
  out
}
