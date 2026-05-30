test_that("sample_cumulants agree with EstemPMM skewness/kurtosis", {
  skip_if_not_installed("EstemPMM")
  set.seed(10)
  x <- rexp(2000) - 1
  sc <- sample_cumulants(x)
  expect_equal(sc$gamma3, EstemPMM::pmm_skewness(x), tolerance = 1e-8)
  expect_equal(sc$gamma4, EstemPMM::pmm_kurtosis(x, excess = TRUE), tolerance = 1e-8)
})

test_that("g2/g3 efficiency factors match the canonical formulas", {
  # asymmetric: g2 = 1 - gamma3^2/(2+gamma4)
  expect_equal(gmdhpmm:::.g2_factor(1.0, 0.5), 1 - 1 / 2.5)
  # platykurtic: g3 = 1 - gamma4^2/(6+9 gamma4+gamma6), positive denominator
  expect_equal(gmdhpmm:::.g3_factor(-1.0, 5.0), 1 - 1 / (6 - 9 + 5))
  # non-positive denominator is guarded to g3 = 1
  expect_equal(gmdhpmm:::.g3_factor(-1.0, 2.0), 1)
  # guarded against non-positive denominator
  expect_equal(gmdhpmm:::.g2_factor(1.0, -3.0), 1)
})

test_that("bootstrap diagnostics are stable and return robust SEs", {
  set.seed(11)
  eps <- rexp(400) - 1
  d <- bootstrap_cumulant_diag(eps, B = 300, seed = 11)
  expect_true(d$gamma3 > 1.0)            # exponential is strongly right-skewed
  expect_true(d$se_gamma3 > 0)
  expect_true(d$se_gamma4 > 0)
  expect_equal(d$B, 300)
  expect_true(d$g2 < 1)                  # PMM2 gain expected
  expect_true(is.finite(d$g2) && is.finite(d$g3))
})

test_that("bootstrap falls back to point estimate on tiny samples", {
  d <- bootstrap_cumulant_diag(rnorm(5), B = 500)
  expect_equal(d$B, 0L)
  expect_true(is.finite(d$se_gamma3))
})
