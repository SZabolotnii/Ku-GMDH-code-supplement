# Regression tests for the GR-transition-fallback hardening (Codex adversarial-review
# finding): missing values in the requested transition channel must NOT silently or
# unsafely change the validation split. Data-independent (builds a synthetic CSV).

.make_volve_csv <- function(n = 600L, gr_na = 0L, all_na = FALSE) {
  set.seed(11)
  dept <- sort(runif(n, 1300, 2500))
  gr <- 40 + 30 * (dept > median(dept)) + rnorm(n)           # a clean step in GR
  if (all_na) gr <- rep(NA_real_, n) else if (gr_na > 0) gr[sample.int(n, gr_na)] <- NA_real_
  d <- data.frame(
    stick_rt = runif(n, 0, 50), rpm = rep(60, n), swob = rnorm(n),
    tqa = 8 + rnorm(n), rop5 = pmax(rexp(n), 1e-3), arc_gr_rt = gr, dept = dept,
    bvel = rnorm(n), source_file = "WL_RAW_BHPR-GR-MECH_TIME_MWD_1.LAS",
    stringsAsFactors = FALSE)
  p <- tempfile(fileext = ".csv"); utils::write.csv(d, p, row.names = FALSE); p
}

test_that("GR transition is used (with provenance) when coverage is adequate, even with sparse NAs", {
  for (na in c(0L, 1L, 30L)) {                                # 30/600 = 5% missing, still >= 80% coverage
    res <- load_volve_drilling(.make_volve_csv(gr_na = na), segment = "gr", k = 100L)
    expect_true(is.finite(res$transition) && length(res$transition) == 1L)
    expect_true(res$transition >= 1L && res$transition <= res$n)
    expect_identical(res$transition_source, "arc_gr_rt")      # did NOT silently fall back
  }
})

test_that("an all-missing GR channel falls back to a finite-safe proxy (no crash, recorded source)", {
  res <- load_volve_drilling(.make_volve_csv(all_na = TRUE), segment = "gr", k = 100L)
  expect_true(is.finite(res$transition) && length(res$transition) == 1L)
  expect_match(res$transition_source, "torque_rop_proxy")     # explicit fallback, not silent GR
})

test_that("the torque/ROP proxy is finite-safe under missing inputs (no NA transition)", {
  p <- .make_volve_csv(all_na = TRUE)
  d <- utils::read.csv(p); d$tqa[c(5, 50, 200)] <- NA; d$rop5[c(10, 60)] <- NA
  p2 <- tempfile(fileext = ".csv"); utils::write.csv(d, p2, row.names = FALSE)
  res <- load_volve_drilling(p2, segment = "gr", k = 100L)
  expect_true(is.finite(res$transition) && length(res$transition) == 1L)
})
