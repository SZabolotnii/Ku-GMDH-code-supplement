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
#   score      U(theta) = Z' [ w (h1 r_c + h2 r_c^2) ] = 0   (r_c centred)
# where (h1, h2) solve the degree-2 Kunchenko normal system F h = b for the
# centred basis {r_c, r_c^2 - c2^w}:
#   F = [[ c2^w, c3^w ], [ c3^w, c4^w + 2 (c2^w)^2 ]],   b = (1, 0)^T
# (b2 = E[d/dtheta (r_c^2)] = E[-2 r_c] = 0 at the weak centre). We do NOT carry
# h1,h2 separately: factor a1 out and run IRLS with weight W_i = w_i (1 + a2 r_ic),
# so the only coefficient needed is the RATIO
#   a2 := h2 / h1 = -c3^w / (2 (c2^w)^2 + c4^w)                         (see l.~50)
# i.e. the regression slope of r_c^2 on r_c. NB this ratio has NO -(g3^w)^2 term:
# the Schur/coupling term that appears in the *stand-alone* coefficient h2 cancels
# in h2/h1 precisely because b2 = 0 (verified numerically against solve(F,b)).
# Do not "restore" a -(g3^w)^2 factor here -- that would be the stand-alone h2,
# not the IRLS ratio this score is factored on. Solved by Gauss-Newton warm-
# started from LSE; weights re-evaluated each step (a redescending, window-
# parametrized estimator). Degree-1 (a2 = 0) reduces to a Welsch M-estimator;
# the a2 term is the skewness-exploiting part.

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

# PATP-3 inner estimator: FRACTIONAL-power signed-parity PMM3 for symmetric heavy
# tails (Kunchenko PATP basis, signed-parity Form B; see the kunchenko-patp skill).
# Where the windowed integer cubic (fit_kg2_weak_pmm3) needs a windowed 6th moment
# m6^w, the fractional route uses only nu_{2p} with 2p < 4 -- finite/stable when the
# residual is symmetric AND heavy (gamma3 ~ 0, large gamma4), exactly the regime a
# well-specified GMDH cascade leaves. Exponent p_3(alpha) = 1/3 + (8/3) alpha^2
# (p_3(0)=1/3 fractal, p_3(1)=3 signed-cube). alpha* is chosen per fit by minimizing
# the closed-form variance-reduction g2(p) on the warm-start residuals; the estimating
# equation sum_i [h1 r_i + h2 sign(r_i)|r_i|^p] z_i = 0 is solved by factored IRLS
# (weight h1 + h2|r_c|^{p-1}, |r_c| floored), h = F2^{-1} (1, p nu_{p-1}).
# NB: provisional extension of the verified S=2/i=2 PATP apparatus (i=3 closed form
# and regression solver are ours; the skill freezes i=2 / location).

.p3_exp <- function(alpha) 1/3 + (8/3) * alpha^2

.g2_from_p <- function(p, c2, r) {
  if (p - 1 <= -1) return(NA_real_)
  a <- abs(r)
  nu_pm1 <- mean(a^(p - 1)); nu_pp1 <- mean(a^(p + 1)); nu_2p <- mean(a^(2 * p)); sig_p <- mean(sign(r) * a^p)
  if (any(!is.finite(c(nu_pm1, nu_pp1, nu_2p, sig_p)))) return(NA_real_)
  F22 <- nu_2p - sig_p^2; num <- c2 * F22 - nu_pp1^2
  den <- c2 * (F22 - 2 * p * nu_pp1 * nu_pm1 + p^2 * c2 * nu_pm1^2)
  if (abs(den) < 1e-12) return(NA_real_); num / den
}

.patp3_alpha_star <- function(r) {
  r <- r - stats::median(r); c2 <- stats::var(r)
  grid <- setdiff(seq(0.05, 0.95, by = 0.05), c(0.45, 0.50, 0.55))
  g <- vapply(grid, function(al) .g2_from_p(.p3_exp(al), c2, r), numeric(1))
  ok <- is.finite(g) & g > 0
  if (!any(ok)) return(0.2)
  grid[ok][which.min(g[ok])]
}

#' Fit a KG-2 partial model by fractional-power PATP-3 (symmetric heavy tails)
#'
#' @param d data.frame with KG-2 features and response y.
#' @param alpha PATP exponent control in (0,1); if NULL (default) chosen per fit by
#'   minimizing the closed-form g2 on the warm-start residuals.
#' @param max_iter,tol IRLS controls.
#' @return length-6 coefficient vector in canonical order.
#' @export
fit_kg2_patp3 <- function(d, alpha = NULL, max_iter = 50L, tol = 1e-7) {
  Z <- stats::model.matrix(.kg2_formula, data = d)
  y <- d$y
  theta <- .fit_kg2_ridge(d, lambda = 1e-8)        # LSE warm start
  if (is.null(alpha)) {
    r0 <- as.numeric(y - Z %*% theta)
    alpha <- .patp3_alpha_star(r0)
  }
  if (abs(alpha - 0.5) < 0.04) return(stats::setNames(as.numeric(theta), .kg2_coef_order))
  p <- .p3_exp(alpha)
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - Z %*% theta)
    s <- stats::mad(r, constant = 1.4826)
    if (!is.finite(s) || s <= 1e-12) s <- stats::sd(r); if (!is.finite(s) || s <= 1e-12) s <- 1
    a <- pmax(abs(r), 0.1 * s)
    nu_pm1 <- mean(a^(p - 1)); nu_pp1 <- mean(a^(p + 1)); nu_2p <- mean(a^(2 * p)); sig_p <- mean(sign(r) * a^p); c2 <- mean(r^2)
    F2 <- matrix(c(c2, nu_pp1, nu_pp1, nu_2p - sig_p^2), 2)
    h <- tryCatch(solve(F2, c(1, p * nu_pm1)), error = function(e) NULL)
    if (is.null(h) || !all(is.finite(h))) break
    W <- h[1] + h[2] * a^(p - 1)
    W[!is.finite(W) | W < 0] <- 0
    if (sum(W) <= 0) break
    theta_new <- .solve_weighted_ridge(Z, y, W, lambda = 1e-8)
    if (any(!is.finite(theta_new))) break
    if (sqrt(sum((theta_new - theta)^2)) / max(1, sqrt(sum(theta^2))) < tol) { theta <- theta_new; break }
    theta <- theta_new
  }
  out <- stats::setNames(as.numeric(theta), .kg2_coef_order)
  out[!is.finite(out)] <- 0
  out
}

# Weak-moment PMM3 inner estimator (symmetric leptokurtic heavy tails). Classical
# PMM3 uses the raw 6th cumulant kappa = (m6 - 3 m4 m2)/(m4 - 3 m2^2), which under
# heavy tails is destabilized (raw m6 may not exist). Weak moments make m6^w finite
# even under Student-t/Cauchy provided the window decays fast (Gaussian / Schwartz),
# so the windowed cubic correction is well defined.
#
# Construction (factored IRLS, mirrors fit_kg2_weak_pmm2):
#   window     w_i = exp(-0.5 ((r_i - mu_w)/(sigma_mult * s))^2),  s = MAD/0.6745
#   weak moms  m2^w, m4^w, m6^w (central, around the weak mean), kappa^w as above
#   weight     W_i = w_i (kappa^w - rc_i^2) / |kappa^w|   (cubic kurtosis weight)
#   score      U(theta) = Z' [ w (kappa^w rc - rc^3) ] = 0
# Weights are redescending (clipped >= 0: residuals beyond |rc| = sqrt(kappa^w) are
# dropped). When the windowed law is mesokurtic (m4^w ~ 3 m2^w^2) kappa^w -> inf and
# it reduces to the degree-1 Welsch M-estimator (W = w). Its practical value is
# bandwidth-insurance: efficiency stays flat at over-wide windows where degree-1
# degrades (see paper-4 run_weak_pmm3_probe.R). NB: this targets symmetric
# LEPTOKURTIC heavy tails; platykurtic residuals (gamma4 < 0) are the classical raw
# PMM3 regime (EstemPMM::lm_pmm3), not this estimator.

#' Fit a KG-2 partial model by weak-moment PMM3
#'
#' @param d data.frame with KG-2 features (b1,b2,b12,b11,b22) and response y.
#' @param sigma_mult window scale as a multiple of the robust residual scale
#'   (default 2.5). Wider windows are tolerated thanks to the cubic correction.
#' @param max_iter,tol IRLS controls.
#' @return length-6 coefficient vector in canonical order.
#' @export
fit_kg2_weak_pmm3 <- function(d, sigma_mult = 2.5, max_iter = 50L, tol = 1e-7) {
  Z <- stats::model.matrix(.kg2_formula, data = d)
  y <- d$y
  theta <- .fit_kg2_ridge(d, lambda = 1e-8)        # LSE warm start
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - Z %*% theta)
    mu <- stats::median(r)
    s <- stats::mad(r, constant = 1.4826)
    if (!is.finite(s) || s <= 1e-12) s <- stats::sd(r)
    if (!is.finite(s) || s <= 1e-12) s <- 1
    sw <- sigma_mult * s
    w <- exp(-0.5 * ((r - mu) / sw)^2)             # Gaussian (Schwartz) window -> finite m6^w
    sw_sum <- sum(w)
    if (!is.finite(sw_sum) || sw_sum <= 0) break
    rc <- r - sum(w * r) / sw_sum                  # centre on weak mean
    m2 <- sum(w * rc^2) / sw_sum
    m4 <- sum(w * rc^4) / sw_sum
    m6 <- sum(w * rc^6) / sw_sum
    dk <- m4 - 3 * m2^2
    if (!is.finite(dk) || abs(dk) < 1e-3 * max(m2^2, 1e-12)) {
      W <- w                                       # mesokurtic -> degree-1 (Welsch)
    } else {
      kappa <- (m6 - 3 * m4 * m2) / dk
      W <- w * (kappa - rc^2) / max(abs(kappa), 1e-8)   # cubic kurtosis weight in (-Inf, 1]
    }
    W[!is.finite(W) | W < 0] <- 0                  # redescending: never negative weight
    if (sum(W) <= 0) W <- w                        # safety: all clipped -> fall back to window
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
