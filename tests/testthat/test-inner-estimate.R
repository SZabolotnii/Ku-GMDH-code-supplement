test_that("inner_estimate returns 6 canonical coefficients", {
  set.seed(30)
  v1 <- rnorm(300); v2 <- rnorm(300)
  y <- 1 + 2 * v1 - v2 + 0.5 * v1 * v2 + rnorm(300)
  est <- inner_estimate(v1, v2, y, gmdh_pmm_control(B = 200, boot_seed = 1))
  expect_length(est$theta, 6)
  expect_named(est$theta, gmdhpmm:::.kg2_coef_order)
  expect_true(est$method %in% c("LSE", "PMM2", "PMM3"))
})

test_that("Gaussian noise dispatches to LSE (safe fallback, H3)", {
  set.seed(31)
  v1 <- rnorm(400); v2 <- rnorm(400)
  y <- 1 + 2 * v1 - v2 + rnorm(400)
  est <- inner_estimate(v1, v2, y, gmdh_pmm_control(B = 200, boot_seed = 2))
  expect_equal(est$method, "LSE")
})

test_that("strongly skewed noise dispatches to PMM2 and refits", {
  set.seed(32)
  v1 <- rnorm(400); v2 <- rnorm(400)
  y <- 1 + 2 * v1 - v2 + (rexp(400) - 1) * 2   # strong right skew
  est <- inner_estimate(v1, v2, y, gmdh_pmm_control(B = 200, boot_seed = 3))
  expect_equal(est$method, "PMM2")
  expect_true(est$converged)
})

test_that("forced ridge-LSE baseline returns canonical coefficients", {
  set.seed(33)
  v1 <- rnorm(120); v2 <- v1 + rnorm(120, sd = 1e-4)
  y <- 1 + 2 * v1 - v2 + rnorm(120, sd = 0.1)
  est <- inner_estimate(v1, v2, y,
                        gmdh_pmm_control(B = 0, force_method = "ridge-LSE",
                                         ridge_lambda = 1e-8))
  expect_equal(est$method, "ridge-LSE")
  expect_length(est$theta, 6)
  expect_named(est$theta, gmdhpmm:::.kg2_coef_order)
  expect_true(all(is.finite(est$theta)))
})

test_that("forced Huber and L1 baselines return canonical coefficients", {
  set.seed(34)
  v1 <- rnorm(160); v2 <- rnorm(160)
  y <- 1 + 2 * v1 - v2 + rnorm(160, sd = 0.2)
  y[seq(1, 160, by = 20)] <- y[seq(1, 160, by = 20)] + 20

  for (method in c("Huber", "L1")) {
    est <- inner_estimate(v1, v2, y,
                          gmdh_pmm_control(B = 0, force_method = method,
                                           max_iter = 40))
    expect_equal(est$method, method)
    expect_length(est$theta, 6)
    expect_named(est$theta, gmdhpmm:::.kg2_coef_order)
    expect_true(all(is.finite(est$theta)))
  }
})

test_that("Huber/L1 reduce leverage from response outliers", {
  set.seed(35)
  v1 <- rnorm(240); v2 <- rnorm(240)
  theta <- c(0.5, 1.5, -1.0, 0.8, 0.6, -0.4)
  y_clean <- kg2_predict(theta, v1, v2) + rnorm(240, sd = 0.2)
  y_bad <- y_clean
  y_bad[seq(1, 240, by = 24)] <- y_bad[seq(1, 240, by = 24)] + 35

  lse <- inner_estimate(v1, v2, y_bad, gmdh_pmm_control(B = 0, force_method = "LSE"))$theta
  huber <- inner_estimate(v1, v2, y_bad,
                          gmdh_pmm_control(B = 0, force_method = "Huber",
                                           max_iter = 60))$theta
  l1 <- inner_estimate(v1, v2, y_bad,
                       gmdh_pmm_control(B = 0, force_method = "L1",
                                        max_iter = 60))$theta

  expect_lt(sum((huber - theta)^2), sum((lse - theta)^2))
  expect_lt(sum((l1 - theta)^2), sum((lse - theta)^2))
})
