# Dispatch sanity, mirroring the EstemPMM dispatch self-test but on the
# bootstrap-stabilized diagnostics (algorithm-spec.md 4.1).

dispatch_on <- function(eps, B = 300, seed = 1) {
  dispatch_method(bootstrap_cumulant_diag(eps, B = B, seed = seed))
}

test_that("Gaussian residuals -> LSE", {
  set.seed(20)
  expect_equal(dispatch_on(rnorm(400)), "LSE")
})

test_that("strongly asymmetric residuals -> PMM2", {
  set.seed(21)
  expect_equal(dispatch_on(rexp(400) - 1), "PMM2")
})

test_that("symmetric platykurtic residuals -> PMM3", {
  set.seed(22)
  # uniform: gamma4 ~ -1.2, symmetric
  eps <- runif(500, -sqrt(3), sqrt(3))
  expect_equal(dispatch_on(eps), "PMM3")
})

test_that("significance gate keeps dispatch safe at moderate sample sizes", {
  # Gaussian residuals must overwhelmingly resolve to LSE once n is large
  # enough for the gate to bite (small-sample kurtosis bias is a known weak
  # spot for the PMM3 path; see algorithm-spec.md 7).
  set.seed(100)
  decisions <- replicate(60, dispatch_method(bootstrap_cumulant_diag(rnorm(120), B = 200)))
  expect_gte(mean(decisions == "LSE"), 0.85)
})
