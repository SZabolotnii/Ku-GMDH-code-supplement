test_that("reserve-aware criterion penalizes under-reserve more than over-reserve", {
  y <- c(100, 120, 200, 600)
  under <- y * 0.8
  over <- y * 1.2
  mse_under <- external_criterion(y, under, type = "MSE")
  mse_over <- external_criterion(y, over, type = "MSE")
  expect_equal(mse_under, mse_over)

  cr_under <- external_criterion(y, under, type = "reserve-aware",
                                 reserve_weight = 2, tail_weight = 2,
                                 overreserve_weight = 0.25)
  cr_over <- external_criterion(y, over, type = "reserve-aware",
                                reserve_weight = 2, tail_weight = 2,
                                overreserve_weight = 0.25)
  expect_gt(cr_under, cr_over)
})

test_that("reserve-aware criterion remains finite on degenerate validation sums", {
  y <- c(-1, 0, 1)
  yhat <- c(0, 0, 0)
  cr <- external_criterion(y, yhat, type = "reserve-aware")
  expect_true(is.finite(cr))
})

test_that("spike-aware criterion upweights under-predicted tail points", {
  y <- c(rep(1, 20), 10, 12)
  under_spike <- y
  over_spike <- y
  under_spike[21:22] <- c(7, 9)
  over_spike[21:22] <- c(13, 15)

  mse_under <- external_criterion(y, under_spike, type = "MSE")
  mse_over <- external_criterion(y, over_spike, type = "MSE")
  expect_equal(mse_under, mse_over)

  cr_under <- external_criterion(y, under_spike, type = "spike-aware",
                                 reserve_weight = 1, tail_weight = 2,
                                 overreserve_weight = 0.25, tail_prob = 0.9)
  cr_over <- external_criterion(y, over_spike, type = "spike-aware",
                                reserve_weight = 1, tail_weight = 2,
                                overreserve_weight = 0.25, tail_prob = 0.9)
  expect_gt(cr_under, cr_over)
})
