
# tests/testthat/test-detrend_temperature.R
#
# Tests for detrend_temperature() and the internal .detrend_single() workhorse.
# Validates trend removal accuracy, multi-column handling, and edge cases.

test_that("loess detrending removes annual trend, preserves diel amplitude", {

  # 90 days hourly: annual trend + diel oscillation
  t_sec <- seq(0, by = 3600, length.out = 90 * 24)
  annual <- 10 * cos(2 * pi * t_sec / (365.25 * 86400))
  diel   <- 3.0 * cos(2 * pi * t_sec / 86400)
  z <- zoo::zoo(10 + annual + diel, order.by = t_sec)

  detrended <- detrend_temperature(z, method = "loess", span = 0.1)

  # Detrended signal should have near-zero mean (trend removed)
  expect_lt(abs(mean(zoo::coredata(detrended), na.rm = TRUE)), 0.5)

  # Diel amplitude should be preserved (~3.0)
  # Check range of detrended signal is close to 2 * diel amplitude
  obs_range <- diff(range(zoo::coredata(detrended), na.rm = TRUE))
  expect_gt(obs_range, 5.0)   # 2 * 3.0 = 6.0, allow some loess bleed
  expect_lt(obs_range, 7.0)
})


test_that("moving average detrending works with default window", {

  # 30 days hourly with linear trend + diel
  t_sec <- seq(0, by = 3600, length.out = 30 * 24)
  trend <- 0.001 * t_sec / 86400  # ~0.03 degC/day
  diel  <- 2.0 * cos(2 * pi * t_sec / 86400)
  z <- zoo::zoo(15 + trend + diel, order.by = t_sec)

  detrended <- detrend_temperature(z, method = "ma", window = 24L)

  # 24-point MA on hourly data averages out the diel cycle,
  # so the trend estimate is the daily mean. Detrended should
  # retain diel oscillation. Check that it oscillates.
  expect_gt(stats::sd(zoo::coredata(detrended), na.rm = TRUE), 1.0)
})


test_that("polynomial detrending removes quadratic trend", {

  t_sec <- seq(0, by = 3600, length.out = 60 * 24)
  # Quadratic trend: peaks at day 30
  trend <- -0.005 * ((t_sec / 86400) - 30)^2 + 20
  diel  <- 1.5 * cos(2 * pi * t_sec / 86400)
  z <- zoo::zoo(trend + diel, order.by = t_sec)

  detrended <- detrend_temperature(z, method = "polynomial", degree = 2L)

  # Quadratic trend should be well removed; residual should be
  # dominated by the diel signal with SD close to diel amplitude / sqrt(2)
  expect_lt(abs(mean(zoo::coredata(detrended), na.rm = TRUE)), 0.5)
  expect_gt(stats::sd(zoo::coredata(detrended), na.rm = TRUE), 0.8)
})


test_that("return_trend = TRUE returns both components", {

  t_sec <- seq(0, by = 3600, length.out = 14 * 24)
  z <- zoo::zoo(sin(2 * pi * t_sec / 86400) + 15, order.by = t_sec)

  result <- detrend_temperature(z, method = "loess", span = 0.3,
                                return_trend = TRUE)

  expect_type(result, "list")
  expect_named(result, c("detrended", "trend"))
  expect_true(zoo::is.zoo(result$detrended))
  expect_true(zoo::is.zoo(result$trend))

  # Trend + detrended should reconstruct original (within floating point)
  reconstructed <- zoo::coredata(result$trend) + zoo::coredata(result$detrended)
  expect_equal(unname(reconstructed), unname(zoo::coredata(z)), tolerance = 1e-10)
})


test_that("multi-column zoo is detrended column-wise", {

  t_sec <- seq(0, by = 3600, length.out = 30 * 24)
  col1 <- 15 + 5 * cos(2 * pi * t_sec / (365.25 * 86400)) +
           2 * cos(2 * pi * t_sec / 86400)
  col2 <- 10 + 3 * cos(2 * pi * t_sec / (365.25 * 86400)) +
           1 * cos(2 * pi * t_sec / 86400)
  z <- zoo::zoo(cbind(shallow = col1, deep = col2), order.by = t_sec)

  detrended <- detrend_temperature(z, method = "loess", span = 0.2)

  expect_true(zoo::is.zoo(detrended))
  expect_equal(ncol(detrended), 2)
  # Both columns should have reduced variance compared to original
  expect_lt(stats::sd(zoo::coredata(detrended)[, 1], na.rm = TRUE),
            stats::sd(zoo::coredata(z)[, 1]))
})


test_that("detrend_temperature rejects non-zoo input", {
  expect_error(detrend_temperature(1:100), "must be a zoo object")
})

