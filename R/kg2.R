# KG-2 partial model: the order-2 Kolmogorov-Gabor polynomial in two variables
#   P(v1, v2; theta) = th0 + th1 v1 + th2 v2 + th3 v1 v2 + th4 v1^2 + th5 v2^2
# (problem-statement.md 3.1). The coefficient vector theta is always stored in
# the canonical order below so the cascade can predict deterministically.

# Canonical coefficient order (matches the lm() design built by kg2_design()).
.kg2_coef_order <- c("(Intercept)", "b1", "b2", "b12", "b11", "b22")

# Formula used for every partial-model fit (LSE and PMM share it).
.kg2_formula <- stats::as.formula(y ~ b1 + b2 + b12 + b11 + b22)

#' KG-2 design features for a variable pair
#'
#' Build the five non-intercept KG-2 features for a pair of input vectors, in
#' the canonical order \code{b1, b2, b12, b11, b22}.
#'
#' @param v1,v2 numeric vectors of equal length.
#' @return a data.frame with columns \code{b1, b2, b12, b11, b22}.
#' @export
kg2_features <- function(v1, v2) {
  data.frame(b1 = v1, b2 = v2, b12 = v1 * v2, b11 = v1^2, b22 = v2^2)
}

#' KG-2 design matrix (features + optional response)
#'
#' @param v1,v2 numeric input vectors.
#' @param y optional response vector; if supplied, added as column \code{y}.
#' @return a data.frame ready for \code{lm()} / \code{EstemPMM::lm_pmm2()}.
#' @export
kg2_design <- function(v1, v2, y = NULL) {
  d <- kg2_features(v1, v2)
  if (!is.null(y)) d$y <- y
  d
}

#' Evaluate a KG-2 partial model
#'
#' @param theta length-6 coefficient vector in canonical order.
#' @param v1,v2 numeric input vectors.
#' @return numeric vector of predictions.
#' @export
kg2_predict <- function(theta, v1, v2) {
  theta[1] + theta[2] * v1 + theta[3] * v2 +
    theta[4] * (v1 * v2) + theta[5] * v1^2 + theta[6] * v2^2
}

# Extract the 6 KG-2 coefficients in canonical order from a fitted model.
# Handles both lm (S3, named coefficients -> reorder by name) and the EstemPMM
# PMM2fit/PMM3fit S4 objects (slot `coefficients`, already in design order but
# unnamed). Coefficients dropped for collinearity (NA) are treated as 0.
.extract_kg2_coef <- function(fit) {
  if (isS4(fit)) {
    out <- as.numeric(methods::slot(fit, "coefficients"))   # canonical design order
    names(out) <- .kg2_coef_order
  } else {
    out <- stats::coef(fit)[.kg2_coef_order]                # lm: reorder by name
  }
  out[is.na(out)] <- 0
  stats::setNames(as.numeric(out), .kg2_coef_order)
}

# Ridge fit for KG-2 baselines. The intercept is left unpenalized so the
# baseline isolates coefficient shrinkage from a mean-shift artifact.
.fit_kg2_ridge <- function(d, lambda = 1e-8) {
  X <- stats::model.matrix(.kg2_formula, data = d)
  y <- d$y
  theta <- .solve_weighted_ridge(X, y, rep(1, length(y)), lambda)
  stats::setNames(theta, .kg2_coef_order)
}

.solve_weighted_ridge <- function(X, y, w, lambda = 1e-8) {
  w <- pmax(as.numeric(w), 0)
  if (!any(w > 0)) w[] <- 1
  w <- w / mean(w[w > 0])
  Xw <- X * sqrt(w)
  yw <- y * sqrt(w)
  penalty <- diag(c(0, rep(lambda, ncol(X) - 1L)), ncol(X))
  tryCatch(
    as.numeric(solve(crossprod(Xw) + penalty, crossprod(Xw, yw))),
    error = function(e) {
      pen_diag <- c(0, rep(lambda, ncol(X) - 1L))
      if (lambda > 0) {
        X_aug <- rbind(Xw, diag(sqrt(pen_diag), ncol(X)))
        y_aug <- c(yw, rep(0, ncol(X)))
        coef <- stats::lm.fit(x = X_aug, y = y_aug)$coefficients
      } else {
        coef <- stats::lm.fit(x = Xw, y = yw)$coefficients
      }
      coef[is.na(coef)] <- 0
      as.numeric(coef)
    }
  )
}

.robust_scale <- function(r) {
  s <- stats::mad(r, constant = 1.4826)
  if (!is.finite(s) || s <= 1e-12) s <- stats::sd(r)
  if (!is.finite(s) || s <= 1e-12) s <- 1
  s
}

# IRLS robust KG-2 baselines. "Huber" is a standard M-estimator; "L1" is LAD
# approximated by IRLS with bounded inverse-residual weights.
.fit_kg2_irls <- function(d, loss = c("Huber", "L1"), huber_k = 1.345,
                          max_iter = 50L, tol = 1e-6, lambda = 1e-8) {
  loss <- match.arg(loss)
  X <- stats::model.matrix(.kg2_formula, data = d)
  y <- d$y
  theta <- .fit_kg2_ridge(d, lambda = lambda)
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - X %*% theta)
    s <- .robust_scale(r)
    u <- abs(r) / s
    if (loss == "Huber") {
      w <- ifelse(u <= huber_k, 1, huber_k / pmax(u, 1e-8))
    } else {
      w <- 1 / pmax(u, 1e-4)
    }
    theta_new <- .solve_weighted_ridge(X, y, w, lambda)
    denom <- max(1, sqrt(sum(theta^2)))
    if (sqrt(sum((theta_new - theta)^2)) / denom < tol) {
      theta <- theta_new
      break
    }
    theta <- theta_new
  }
  stats::setNames(as.numeric(theta), .kg2_coef_order)
}
