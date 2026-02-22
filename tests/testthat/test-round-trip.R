
# tests/testthat/test-round-trip.R
#
# Round-trip validation: forward model -> inverse model -> compare.
# Tests that amplitude/phase extraction methods correctly recover known
# thermal parameters from synthetic temperature time series.

test_that("OLS round-trip recovers parameters to machine precision", {
  results <- test_round_trip(fit_method = "ols")

  # OLS on perfect cosines should recover to floating-point precision
  expect_true(
    all(results$relative_error < 1e-10),
    info = paste("Max relative error:", max(results$relative_error))
  )
})

test_that("NLS round-trip recovers parameters within noise tolerance", {
  results <- test_round_trip(fit_method = "nls")

  # nls injects ~0.01% noise to avoid singular gradient on perfect cosines,

  # so we expect ~1e-3 relative error at worst
  expect_true(
    all(results$relative_error < 1e-2),
    info = paste("Max relative error:", max(results$relative_error))
  )
})

test_that("OLS and NLS agree on amplitude ratios within noise tolerance", {
  ser <- generate_example_series()
  aquifer <- generate_example_aquifer()

  # Pick a representative case: diel, first parameter combo
  ts_data <- ser$diel[[1]]$timeSeries
  signal  <- ser$diel[[1]]$signal

  obs_ols <- thObservedSeries(
    empiricalData = ts_data,
    xVals         = ser$diel[[1]]$xVals,
    aquifer       = aquifer,
    boundaryMean  = signal$boundary$mean,
    period        = signal$boundary$period,
    headGrad      = signal$hydro$headGrad,
    nmin          = 10,
    fit_method    = "ols"
  )

  obs_nls <- thObservedSeries(
    empiricalData = ts_data,
    xVals         = ser$diel[[1]]$xVals,
    aquifer       = aquifer,
    boundaryMean  = signal$boundary$mean,
    period        = signal$boundary$period,
    headGrad      = signal$hydro$headGrad,
    nmin          = 10,
    fit_method    = "nls"
  )

  # Amplitude ratios should agree within the nls noise floor
  expect_equal(obs_ols$ampRatio, obs_nls$ampRatio, tolerance = 1e-3)
})

test_that("All scenario/parameter combinations are tested", {
  results <- test_round_trip()

  # 3 diel combos + 2 annual combos = 5 combos, 3 parameters each = 15 rows
  expect_equal(nrow(results), 15)
  expect_equal(sort(unique(results$scenario)), c("annual", "diel"))
  expect_equal(
    sort(unique(results$parameter)),
    c("advectiveThermVel", "darcyFlux", "diffusivity_effective")
  )
})




# --- FFT method tests ---
# FFT requires records where the target period lands near a Fourier bin.
# For diel (period = 86400s, dt = 3600s): n must be a multiple of 24.

test_that("FFT round-trip recovers parameters at machine precision (bin-aligned)", {
  fft_scenarios <- list(
    diel = list(
      k = 10 / 86400, dispersivity = c(0.001, 0.01, 0.1),
      headGrad = 0.1, mean = 22, amplitude = 6,
      phase = 86400 + 3600, period = 86400,
      xVals = c(0, 0.33, 0.66, 1),
      tVals = seq(0, 10 * 86400 - 3600, 3600)  # 240 points, exact bin
    )
  )
  results <- test_round_trip(scenarios = fft_scenarios, fit_method = "fft")
  expect_true(max(results$relative_error) < 1e-10)
})

test_that("FFT agrees with OLS on amplitude ratios (bin-aligned)", {
  fft_scenarios <- list(
    diel = list(
      k = 10 / 86400, dispersivity = c(0.001, 0.01, 0.1),
      headGrad = 0.1, mean = 22, amplitude = 6,
      phase = 86400 + 3600, period = 86400,
      xVals = c(0, 0.33, 0.66, 1),
      tVals = seq(0, 10 * 86400 - 3600, 3600)
    )
  )
  fwd <- generate_example_series(scenarios = fft_scenarios)
  obs_fft <- generate_example_observed(fwd, fit_method = "fft")
  obs_ols <- generate_example_observed(fwd, fit_method = "ols")

  for (i in seq_along(obs_fft$diel)) {
    fft_ar <- obs_fft$diel[[i]]$ampRatio[1, 2]
    ols_ar <- obs_ols$diel[[i]]$ampRatio[1, 2]
    expect_equal(fft_ar, ols_ar, tolerance = 1e-10)
  }
})

test_that("FFT warns on poor bin alignment", {
  short_scenarios <- list(
    diel = list(
      k = 10 / 86400, dispersivity = 0.01,
      headGrad = 0.1, mean = 22, amplitude = 6,
      phase = 86400 + 3600, period = 86400,
      xVals = c(0, 0.5),
      tVals = seq(0, 2 * 86400, 3600)  # 49 points, bin at 24.5h
    )
  )
  fwd <- generate_example_series(scenarios = short_scenarios)
  expect_warning(
    generate_example_observed(fwd, fit_method = "fft"),
    "nearest FFT bin period"
  )
})

