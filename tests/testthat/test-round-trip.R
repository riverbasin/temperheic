
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

