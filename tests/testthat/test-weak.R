# Tests for the weak-moment machinery (Paper 4): WPMM3, weak cumulant diagnostics,
# and the weak dispatch. The classical path is covered by the other test files.

make_kg2 <- function(n, errfn, seed = 1) {
  set.seed(seed)
  v1 <- rnorm(n); v2 <- rnorm(n)
  mu <- 0.5 + 1.0 * v1 - 0.8 * v2 + 0.4 * v1 * v2 - 0.3 * v1^2 + 0.2 * v2^2
  list(v1 = v1, v2 = v2, d = kg2_design(v1, v2, mu + errfn(n)),
       theta = c(0.5, 1.0, -0.8, 0.4, -0.3, 0.2))
}

test_that("fit_kg2_weak_pmm3 returns finite canonical coefficients", {
  m <- make_kg2(300, function(n) rt(n, 3) / sqrt(3))
  th <- fit_kg2_weak_pmm3(m$d)
  expect_length(th, 6)
  expect_true(all(is.finite(th)))
  expect_named(th, .kg2_coef_order)
  expect_lt(sqrt(sum((th - m$theta)^2)), 0.5)   # recovers the truth reasonably
})

test_that("WPMM3 ~ degree-1 weak under Gaussian (mesokurtic reduction)", {
  m <- make_kg2(400, function(n) rnorm(n))
  th3 <- fit_kg2_weak_pmm3(m$d)
  th1 <- fit_kg2_weak_pmm2(m$d)               # a2->~0 under symmetry, ~ degree-1
  expect_lt(sqrt(sum((th3 - th1)^2)), 0.25)   # close under mesokurtic Gaussian
})

test_that("weak cumulants stay finite under heavy tails where raw blow up", {
  set.seed(7)
  eps_cauchy <- rcauchy(500)
  wc <- weak_sample_cumulants(eps_cauchy)
  expect_true(is.finite(wc$gamma3) && is.finite(wc$gamma4) && is.finite(wc$gamma6))
  expect_true(abs(wc$gamma4) < 1e3)           # windowed kurtosis bounded
})

test_that("tail_ratio ~ 1 for Gaussian and large for heavy tails", {
  set.seed(11)
  wg <- weak_cumulant_diag(rnorm(2000), B = 0)
  wh <- weak_cumulant_diag(rt(2000, df = 2), B = 0)
  expect_lt(abs(wg$tail_ratio - 1), 0.25)
  expect_gt(wh$tail_ratio, wg$tail_ratio)
})

test_that("dispatch_method_weak routes by regime", {
  set.seed(3)
  asym_heavy <- { e <- rnorm(4000); h <- runif(4000) < 0.12; e[h] <- -abs(rt(sum(h), 2)) * 2.5; e }
  sym_heavy  <- rt(4000, df = 2)
  light      <- rnorm(4000)
  dw_a <- weak_cumulant_diag(asym_heavy, B = 0); dr_a <- bootstrap_cumulant_diag(asym_heavy, B = 0)
  dw_s <- weak_cumulant_diag(sym_heavy, B = 0);  dr_s <- bootstrap_cumulant_diag(sym_heavy, B = 0)
  dw_l <- weak_cumulant_diag(light, B = 0);      dr_l <- bootstrap_cumulant_diag(light, B = 0)
  expect_equal(dispatch_method_weak(dw_a, dr_a), "WPMM2")
  expect_equal(dispatch_method_weak(dw_s, dr_s), "WPMM3")
  expect_equal(dispatch_method_weak(dw_l, dr_l), "LSE")
})

test_that("inner_estimate supports WPMM3 and auto-weak without breaking classical auto", {
  m <- make_kg2(300, function(n) { e <- rnorm(n); h <- runif(n) < 0.12; e[h] <- -abs(rt(sum(h), 2)) * 2.5; e })
  e_w3 <- inner_estimate(m$v1, m$v2, m$d$y, gmdh_pmm_control(B = 0, force_method = "WPMM3"))
  expect_true(all(is.finite(e_w3$theta)))
  e_aw <- inner_estimate(m$v1, m$v2, m$d$y, gmdh_pmm_control(B = 100, force_method = "auto-weak"))
  expect_true(e_aw$method %in% c("WPMM2", "WPMM3", "LSE", "PMM2", "PMM3"))
  expect_false(is.null(e_aw$diag_weak))
  e_auto <- inner_estimate(m$v1, m$v2, m$d$y, gmdh_pmm_control(B = 100))
  expect_true(e_auto$method %in% c("LSE", "PMM2", "PMM3"))   # classical path unaffected
})

test_that("auto-valgate (contiguous k-fold inner CV) runs, honors valgate_folds/candidates, selects a valid method", {
  expect_equal(gmdh_pmm_control(valgate_folds = 3L)$valgate_folds, 3L)   # param propagates
  m <- make_kg2(400, function(n) { e <- rnorm(n); h <- runif(n) < 0.12; e[h] <- -abs(rt(sum(h), 2)) * 2.5; e })
  e_vg <- inner_estimate(m$v1, m$v2, m$d$y,
                         gmdh_pmm_control(B = 0, force_method = "auto-valgate", valgate_folds = 4L))
  expect_true(all(is.finite(e_vg$theta)))
  expect_true(e_vg$method %in% c("LSE", "Huber", "L1", "WPMM2", "WPMM3"))   # picked from the candidate pool
  # a restricted candidate set must be honored by the validation gate
  e_vg2 <- inner_estimate(m$v1, m$v2, m$d$y,
                          gmdh_pmm_control(B = 0, force_method = "auto-valgate",
                                           valgate_candidates = c("LSE", "WPMM2"), valgate_folds = 3L))
  expect_true(e_vg2$method %in% c("LSE", "WPMM2"))
  # the inner-fold ablation: default stays blocked, "random" runs and leaves the RNG untouched
  expect_equal(gmdh_pmm_control()$valgate_inner, "blocked")
  expect_error(gmdh_pmm_control(valgate_inner = "shuffled"))
  set.seed(7); r0 <- runif(1); set.seed(7)
  e_vr <- inner_estimate(m$v1, m$v2, m$d$y,
                         gmdh_pmm_control(B = 0, force_method = "auto-valgate", valgate_inner = "random"))
  expect_equal(runif(1), r0)
  expect_true(e_vr$method %in% c("LSE", "Huber", "L1", "WPMM2", "WPMM3"))
})
