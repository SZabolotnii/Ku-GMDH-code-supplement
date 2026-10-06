test_that("FORGE standard loader removes units row and returns canonical frame", {
  tmp <- tempfile(fileext = ".csv")
  header <- paste(
    "Date", "Bit Depth", "Depth Hole Total Vertical Depth",
    "Rate of Penetration (Depth/Hour)", "Weight on Bit",
    "Rotary Revolutions per Minute", "Differential Pressure",
    "Flow In", "Top Drive Torque", sep = ","
  )
  units <- paste("date/time", "feet (ft)", "feet (ft)",
                 "feet per hour (ft/hr)", "kilopounds (klbs)",
                 "revolutions per minute (RPM)", "psi",
                 "gallons per minute (GPM)", "kilofoot-pounds (kft-lbs)",
                 sep = ",")
  rows <- vapply(seq_len(80), function(i) {
    paste(
      sprintf("2021-02-08 00:%02d:00", i %% 60),
      100 + i,
      110 + i,
      20 + 0.3 * i,
      5 + 0.05 * i,
      60 + i %% 5,
      100 + i,
      400 + i,
      3 + 0.02 * i,
      sep = ","
    )
  }, character(1))
  writeLines(c(header, units, rows), tmp)

  D <- load_forge_standard_drilling(path = tmp, k = 5L)
  expect_equal(D$n, 80)
  expect_equal(D$target, "log_rop")
  expect_true(all(c("wob", "rpm", "dept", "tvd") %in% D$predictors))
  expect_true(all(is.finite(D$X)))
  expect_true(all(is.finite(D$y)))
  expect_true(D$transition >= 1L && D$transition <= D$n)
})

test_that("drilling data manifest records missing and present catalog entries", {
  root <- tempdir()
  dir.create(file.path(root, "forge"), showWarnings = FALSE)
  path <- file.path(root, "forge", "Well_58-32_processed_pason_log.csv")
  writeLines(c("x,y", "1,2", "3,4"), path)

  m <- drilling_data_manifest(root = root, count_rows = TRUE)
  expect_true("forge58_processed" %in% m$dataset_id)
  hit <- m[m$dataset_id == "forge58_processed", ]
  expect_true(hit$exists)
  expect_equal(hit$line_count, 3L)
  expect_true(any(!m$exists))
})
