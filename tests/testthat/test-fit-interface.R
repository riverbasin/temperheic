
# tests/testthat/test-fit-interface.R
#
# Verifies that all three fitting methods (fit_ols, fit_fft, fitCosine)
# produce identical output structures, so thObservedSeries can use them
# as drop-in replacements.

test_that("fit_ols and fit_fft return identical structure", {

  # 10 days hourly, bin-aligned for FFT
  t_sec <- seq(0, by = 3600, length.out = 10 * 24)
  vals <- cbind(
    x0   = 6 * cos(2 * pi * t_sec / 86400) + 22,
    x0.5 = 4 * cos(2 * pi * t_sec / 86400 - 0.3) + 22
  )
  z <- zoo::zoo(vals, order.by = t_sec)

  ols <- fit_ols(z, boundaryMean = 22, periodInSeconds = 86400,
                 optimizeRange = c(-1/8, 7/8), nmin = 10,
                 empiricalDataPeriods = c(1, 1))

  fft_result <- fit_fft(z, boundaryMean = 22, periodInSeconds = 86400,
                        optimizeRange = c(-1/8, 7/8), nmin = 10,
                        empiricalDataPeriods = c(1, 1))

  # Same list components
  expect_named(ols, c("deltaPhaseRadians", "ampRatio"))
  expect_named(fft_result, c("deltaPhaseRadians", "ampRatio"))

  # Same attributes
  expect_true(!is.null(attr(ols, "amplitudes")))
  expect_true(!is.null(attr(ols, "phases")))
  expect_true(!is.null(attr(fft_result, "amplitudes")))
  expect_true(!is.null(attr(fft_result, "phases")))

  # Same lengths (n^2 for 2 sensors = 4)
  expect_length(ols$deltaPhaseRadians, 4)
  expect_length(fft_result$deltaPhaseRadians, 4)

  # Values should agree on bin-aligned data
  expect_equal(attr(ols, "amplitudes"), attr(fft_result, "amplitudes"),
               tolerance = 1e-6)
})


test_that("fitCosine returns same structure as fit_ols", {

  # fitCosine (NLS) needs non-trivial residuals to converge. Add realistic
  # noise to avoid the singular-gradient failure on perfect cosines.
  # This mirrors what fitCosine itself does internally, but we add more
  # noise here to ensure convergence in all test environments.
  set.seed(42)
  t_sec <- seq(0, by = 3600, length.out = 10 * 24)
  noise1 <- rnorm(length(t_sec), sd = 0.1)
  noise2 <- rnorm(length(t_sec), sd = 0.1)
  vals <- cbind(
    x0   = 6 * cos(2 * pi * t_sec / 86400) + 22 + noise1,
    x0.5 = 4 * cos(2 * pi * t_sec / 86400 - 0.3) + 22 + noise2
  )
  z <- zoo::zoo(vals, order.by = t_sec)

  ols <- fit_ols(z, boundaryMean = 22, periodInSeconds = 86400,
                 optimizeRange = c(-1/8, 7/8), nmin = 10,
                 empiricalDataPeriods = c(1, 1))

  nls_result <- fitCosine(z, boundaryMean = 22, periodInSeconds = 86400,
                          optimizeRange = c(-1/8, 7/8), nmin = 10,
                          empiricalDataPeriods = c(1, 1))

  # Identical structure
  expect_named(nls_result, names(ols))
  expect_length(nls_result$ampRatio, length(ols$ampRatio))
  expect_true(!is.null(attr(nls_result, "amplitudes")))
  expect_true(!is.null(attr(nls_result, "phases")))

  # Values should agree within NLS noise tolerance
  expect_equal(attr(ols, "amplitudes"), attr(nls_result, "amplitudes"),
               tolerance = 0.05)
})


test_that("fit methods handle nmin threshold consistently", {

  # Very short record: only 5 observations, below nmin = 10
  t_sec <- seq(0, by = 3600, length.out = 5)
  vals <- cbind(x0 = sin(2 * pi * t_sec / 86400) + 15)
  z <- zoo::zoo(vals, order.by = t_sec)

  ols <- fit_ols(z, boundaryMean = 15, periodInSeconds = 86400,
                 optimizeRange = c(-1/8, 7/8), nmin = 10,
                 empiricalDataPeriods = 1)

  # FFT on 5 points will warn about bin mismatch — expected and harmless
  fft_r <- suppressWarnings(
    fit_fft(z, boundaryMean = 15, periodInSeconds = 86400,
            optimizeRange = c(-1/8, 7/8), nmin = 10,
            empiricalDataPeriods = 1)
  )

  # Both should return NA amplitudes when n < nmin
  expect_true(is.na(attr(ols, "amplitudes")[[1]]))
  expect_true(is.na(attr(fft_r, "amplitudes")[[1]]))
})

