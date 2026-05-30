test_that("kg2_predict matches the explicit KG-2 polynomial", {
  set.seed(1)
  v1 <- rnorm(20); v2 <- rnorm(20)
  theta <- c(0.5, 1, -2, 0.3, 0.7, -0.4)
  manual <- theta[1] + theta[2] * v1 + theta[3] * v2 +
    theta[4] * v1 * v2 + theta[5] * v1^2 + theta[6] * v2^2
  expect_equal(kg2_predict(theta, v1, v2), manual)
})

test_that("LSE recovers known KG-2 coefficients in canonical order", {
  set.seed(2)
  v1 <- rnorm(500); v2 <- rnorm(500)
  theta <- c(0.5, 1, -2, 0.3, 0.7, -0.4)
  y <- kg2_predict(theta, v1, v2) + rnorm(500, sd = 1e-3)
  est <- gmdhpmm:::.extract_kg2_coef(stats::lm(gmdhpmm:::.kg2_formula,
                                               data = kg2_design(v1, v2, y)))
  expect_named(est, gmdhpmm:::.kg2_coef_order)
  expect_equal(unname(est), theta, tolerance = 1e-2)
})

test_that("kg2_design has the canonical feature columns", {
  d <- kg2_design(1:3, 4:6, y = 7:9)
  expect_equal(names(d), c("b1", "b2", "b12", "b11", "b22", "y"))
  expect_equal(d$b12, (1:3) * (4:6))
  expect_equal(d$b11, (1:3)^2)
})

test_that("weighted ridge fallback handles rank-deficient systems", {
  X <- cbind(1, x = 1:6, x_dup = 1:6)
  y <- 1 + 2 * (1:6)
  theta <- gmdhpmm:::.solve_weighted_ridge(X, y, w = rep(1, 6), lambda = 0)
  expect_length(theta, 3)
  expect_true(all(is.finite(theta)))
  expect_equal(drop(X %*% theta), y, tolerance = 1e-8)
})
