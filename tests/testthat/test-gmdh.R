# End-to-end tournament behaviour.

make_data <- function(n, d = 4, noise) {
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  # nonlinear target on a couple of inputs
  y <- 1 + 1.5 * X[, 1] - X[, 2] + 0.8 * X[, 1] * X[, 2] + 0.5 * X[, 3]^2 + noise
  list(X = X, y = y)
}

test_that("gmdh_pmm fits, predicts, and reports a finite criterion", {
  set.seed(40)
  dat <- make_data(400, noise = rnorm(400))
  fit <- gmdh_pmm(dat$X, dat$y,
                  gmdh_pmm_control(B = 100, L_max = 4, F = 6, seed = 40))
  expect_s3_class(fit, "gmdh_pmm")
  expect_true(is.finite(fit$best$cr))

  test <- make_data(200, noise = rnorm(200))
  yhat <- predict(fit, test$X)
  expect_length(yhat, 200)
  expect_true(all(is.finite(yhat)))
  # beats the intercept-only baseline out of sample
  expect_lt(mean((test$y - yhat)^2), stats::var(test$y))
})

test_that("under Gaussian noise the dispatch is overwhelmingly LSE (H3)", {
  # Safe fallback: layer-0 residuals are Gaussian so dispatch picks LSE; deeper
  # layers can drift slightly, so we assert a dominant LSE share rather than all.
  set.seed(41)
  dat <- make_data(500, noise = rnorm(500))
  fit <- gmdh_pmm(dat$X, dat$y,
                  gmdh_pmm_control(B = 150, L_max = 3, F = 6, seed = 41))
  methods <- vapply(Filter(function(z) !is.null(z) && z$type == "model", fit$nodes),
                    function(z) z$method, character(1))
  expect_gte(mean(methods == "LSE"), 0.85)
})

test_that("per-layer best_id supports depth-curve prediction", {
  set.seed(43)
  dat <- make_data(400, noise = rnorm(400))
  fit <- gmdh_pmm(dat$X, dat$y,
                  gmdh_pmm_control(B = 50, L_max = 4, F = 6, epsilon = -Inf, seed = 43))
  expect_gte(length(fit$layers), 2)
  for (L in seq_along(fit$layers)) {
    id <- fit$layers[[L]]$best_id
    expect_false(is.null(id))
    yhat <- predict(fit, dat$X, node = id)
    expect_length(yhat, nrow(dat$X))
    expect_true(all(is.finite(yhat)))
  }
})

test_that("fixed train/validation partitions are respected", {
  set.seed(44)
  dat <- make_data(240, noise = rnorm(240))
  fit <- gmdh_pmm(dat$X, dat$y,
                  gmdh_pmm_control(B = 20, L_max = 2, F = 4,
                                   train_index = 1:160,
                                   val_index = 161:240))
  expect_equal(fit$n_train, 160)
  expect_equal(fit$n_val, 80)
  expect_true(is.finite(fit$best$cr))

  expect_error(
    gmdh_pmm(dat$X, dat$y,
             gmdh_pmm_control(B = 0, train_index = 1:20, val_index = 20:40)),
    "must not overlap"
  )
})

test_that("predict errors on wrong column count", {
  set.seed(42)
  dat <- make_data(200, noise = rnorm(200))
  fit <- gmdh_pmm(dat$X, dat$y, gmdh_pmm_control(B = 50, L_max = 2, seed = 42))
  expect_error(predict(fit, dat$X[, 1:2]), "columns")
})
