# Residual cumulant diagnostics with bootstrap stabilization.
#
# Naive sample cumulants on the small partitions inside a GMDH tournament
# (50-500 points) are very noisy, so the dispatch can flip on sampling noise.
# bootstrap_cumulant_diag() stabilizes them with B resamples and a robust
# (median / MAD) summary -- this is contribution C2 of Paper 1
# (problem-statement.md 7, algorithm-spec.md 7).

# PMM2 efficiency factor g2 = 1 - gamma3^2 / (2 + gamma4)  (algorithm-spec 5.1).
# Guarded: a non-positive denominator means the formula does not apply -> g2 = 1.
.g2_factor <- function(gamma3, gamma4) {
  denom <- 2 + gamma4
  if (!is.finite(denom) || denom <= 0) return(1)
  1 - gamma3^2 / denom
}

# PMM3 efficiency factor g3 = 1 - gamma4^2 / (6 + 9 gamma4 + gamma6) (6.1).
.g3_factor <- function(gamma4, gamma6) {
  denom <- 6 + 9 * gamma4 + gamma6
  if (!is.finite(denom) || denom <= 0) return(1)
  1 - gamma4^2 / denom
}

#' Sample standardized cumulants of a vector
#'
#' Central moments and standardized cumulants up to order 6, using the same
#' conventions as \code{EstemPMM} and the dispatch reference: \code{gamma3}
#' (skewness), \code{gamma4} (excess kurtosis), \code{gamma6} (standardized
#' 6th cumulant \eqn{m6/m2^3 - 15 m4/m2^2 + 30}).
#'
#' @param eps numeric vector (e.g. residuals).
#' @return named list \code{m2, m3, m4, m6, gamma3, gamma4, gamma6}.
#' @export
sample_cumulants <- function(eps) {
  r <- eps - mean(eps)
  m2 <- mean(r^2); m3 <- mean(r^3); m4 <- mean(r^4); m6 <- mean(r^6)
  list(
    m2 = m2, m3 = m3, m4 = m4, m6 = m6,
    gamma3 = m3 / m2^(3 / 2),
    gamma4 = m4 / m2^2 - 3,
    gamma6 = m6 / m2^3 - 15 * (m4 / m2^2) + 30
  )
}

#' Bootstrap-stabilized cumulant diagnostics
#'
#' Resamples \code{eps} with replacement \code{B} times, summarizing the
#' standardized cumulants robustly (median point estimate, MAD-based standard
#' error). Falls back to the point estimate with parametric Gaussian SEs when
#' the sample is too small to bootstrap (\code{n < 8}) or \code{B <= 0}.
#'
#' @param eps numeric residual vector.
#' @param B number of bootstrap resamples (default 500).
#' @param robust if TRUE (default) use median / MAD; else mean / sd.
#' @param seed optional RNG seed for reproducibility.
#' @return named list with \code{gamma3, gamma4, gamma6}, robust SEs
#'   \code{se_gamma3, se_gamma4}, efficiency factors \code{g2, g3}, the sample
#'   size \code{n}, the effective \code{B}, and the raw point estimate
#'   \code{point} (used by the PMM-loss criterion).
#' @export
bootstrap_cumulant_diag <- function(eps, B = 500, robust = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(eps)
  pt <- sample_cumulants(eps)

  if (B <= 0 || n < 8) {
    g3hat <- pt$gamma3; g4hat <- pt$gamma4; g6hat <- pt$gamma6
    se3 <- sqrt(6 / n); se4 <- sqrt(24 / n)         # Gaussian asymptotic SEs
    Beff <- 0L
  } else {
    G3 <- numeric(B); G4 <- numeric(B); G6 <- numeric(B)
    for (b in seq_len(B)) {
      cb <- sample_cumulants(eps[sample.int(n, n, replace = TRUE)])
      G3[b] <- cb$gamma3; G4[b] <- cb$gamma4; G6[b] <- cb$gamma6
    }
    if (robust) {
      g3hat <- stats::median(G3); g4hat <- stats::median(G4); g6hat <- stats::median(G6)
      se3 <- stats::mad(G3, constant = 1.4826); se4 <- stats::mad(G4, constant = 1.4826)
    } else {
      g3hat <- mean(G3); g4hat <- mean(G4); g6hat <- mean(G6)
      se3 <- stats::sd(G3); se4 <- stats::sd(G4)
    }
    Beff <- B
  }

  list(
    gamma3 = g3hat, gamma4 = g4hat, gamma6 = g6hat,
    se_gamma3 = se3, se_gamma4 = se4,
    g2 = .g2_factor(g3hat, g4hat),
    g3 = .g3_factor(g4hat, g6hat),
    n = n, B = Beff, point = pt
  )
}

#' Weak (value-domain windowed) standardized cumulants of a vector
#'
#' Replaces raw moment averaging \eqn{E[g(r)]} by a value-domain windowed
#' functional \eqn{E_w[g(r)] = \sum_i w_\sigma(r_i) g(r_i) / \sum_i w_\sigma(r_i)},
#' with a Gaussian window centred on the robust residual location, then forms the
#' standardized cumulants from the windowed central moments. These stay finite and
#' stable under heavy tails. IMPORTANT: these are cumulants of the windowed
#' functional, NOT true cumulants of the residual distribution.
#'
#' @param eps numeric residual vector.
#' @param sigma_mult window scale as a multiple of the robust scale (default 2.5).
#' @return named list \code{gamma3, gamma4, gamma6} (weak standardized cumulants),
#'   and the weak central moments \code{c2, c3, c4, c6}.
#' @export
weak_sample_cumulants <- function(eps, sigma_mult = 2.5) {
  mu <- stats::median(eps)
  s <- stats::mad(eps, constant = 1.4826)
  if (!is.finite(s) || s <= 1e-12) s <- stats::sd(eps)
  if (!is.finite(s) || s <= 1e-12) s <- 1
  w <- exp(-0.5 * ((eps - mu) / (sigma_mult * s))^2)
  sws <- sum(w)
  if (!is.finite(sws) || sws <= 0) {
    return(list(gamma3 = 0, gamma4 = 0, gamma6 = 0, c2 = NA, c3 = NA, c4 = NA, c6 = NA))
  }
  mw <- sum(w * eps) / sws
  rc <- eps - mw
  c2 <- sum(w * rc^2) / sws
  c3 <- sum(w * rc^3) / sws
  c4 <- sum(w * rc^4) / sws            # windowed central 4th moment (raw, not excess)
  c6 <- sum(w * rc^6) / sws
  if (!is.finite(c2) || c2 <= 1e-12) {
    return(list(gamma3 = 0, gamma4 = 0, gamma6 = 0, c2 = c2, c3 = c3, c4 = c4, c6 = c6))
  }
  list(
    c2 = c2, c3 = c3, c4 = c4, c6 = c6,
    gamma3 = c3 / c2^(3 / 2),
    gamma4 = c4 / c2^2 - 3,
    gamma6 = c6 / c2^3 - 15 * (c4 / c2^2) + 30
  )
}

#' Bootstrap-stabilized WEAK cumulant diagnostics (for the weak dispatch)
#'
#' Weak-moment analogue of \code{\link{bootstrap_cumulant_diag}}: bootstrap-
#' stabilized windowed (weak) standardized cumulants \eqn{\gamma_3^w, \gamma_4^w,
#' \gamma_6^w}, plus a heavy-tail indicator \code{tail_ratio = sd / (MAD/0.6745)}
#' (=1 for Gaussian, large under heavy tails / contamination) and the raw excess
#' kurtosis \code{gamma4_raw} -- both used by \code{\link{dispatch_method_weak}} to
#' gate weak-polynomial corrections to genuinely heavy-tailed regimes.
#'
#' @param eps numeric residual vector.
#' @param B number of bootstrap resamples (default 500).
#' @param sigma_mult window scale multiple (default 2.5).
#' @param robust if TRUE (default) use median / MAD summary; else mean / sd.
#' @param seed optional RNG seed.
#' @return named list with weak \code{gamma3, gamma4, gamma6}, robust SEs
#'   \code{se_gamma3, se_gamma4}, heavy-tail signals \code{tail_ratio, gamma4_raw},
#'   weak efficiency factors \code{g2, g3}, \code{n}, effective \code{B}, and the
#'   weak point estimate \code{point}.
#' @export
weak_cumulant_diag <- function(eps, B = 500, sigma_mult = 2.5, robust = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(eps)
  pt <- weak_sample_cumulants(eps, sigma_mult = sigma_mult)
  raw <- sample_cumulants(eps)
  sd_eps <- stats::sd(eps)
  mad_scale <- stats::mad(eps, constant = 1.4826)
  tail_ratio <- if (is.finite(mad_scale) && mad_scale > 1e-12) sd_eps / mad_scale else NA_real_

  if (B <= 0 || n < 8) {
    g3hat <- pt$gamma3; g4hat <- pt$gamma4; g6hat <- pt$gamma6
    se3 <- sqrt(6 / n); se4 <- sqrt(24 / n)
    Beff <- 0L
  } else {
    G3 <- numeric(B); G4 <- numeric(B); G6 <- numeric(B)
    for (b in seq_len(B)) {
      cb <- weak_sample_cumulants(eps[sample.int(n, n, replace = TRUE)], sigma_mult = sigma_mult)
      G3[b] <- cb$gamma3; G4[b] <- cb$gamma4; G6[b] <- cb$gamma6
    }
    if (robust) {
      g3hat <- stats::median(G3); g4hat <- stats::median(G4); g6hat <- stats::median(G6)
      se3 <- stats::mad(G3, constant = 1.4826); se4 <- stats::mad(G4, constant = 1.4826)
    } else {
      g3hat <- mean(G3); g4hat <- mean(G4); g6hat <- mean(G6)
      se3 <- stats::sd(G3); se4 <- stats::sd(G4)
    }
    Beff <- B
  }

  list(
    gamma3 = g3hat, gamma4 = g4hat, gamma6 = g6hat,
    se_gamma3 = se3, se_gamma4 = se4,
    tail_ratio = tail_ratio, gamma4_raw = raw$gamma4,
    g2 = .g2_factor(g3hat, g4hat),
    g3 = .g3_factor(g4hat, g6hat),
    n = n, B = Beff, point = pt
  )
}
