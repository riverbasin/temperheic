
# tests/testthat/test-signal-generators.R
#
# Tests for generate_composite_boundary() and S3 constructor validation.

# === generate_composite_boundary =============================================

test_that("generate_composite_boundary returns zoo with correct length", {

  z <- generate_composite_boundary(n_years = 1, dt = 3600)

  expect_true(zoo::is.zoo(z))

  # 1 year at hourly: 365.25 * 24 = 8766 observations
  expect_equal(length(z), floor(365.25 * 24))

  # Index should be in seconds, starting at 0
 expect_equal(as.numeric(zoo::index(z)[1]), 0)
})


test_that("generate_composite_boundary has correct spectral content", {

  # 2 years hourly — long enough for clean FFT
  z <- generate_composite_boundary(n_years = 2, dt = 3600)
  n <- length(z)
  dt <- 3600

  # FFT
  fft_result <- stats::fft(zoo::coredata(z) - mean(zoo::coredata(z)))
  power <- Mod(fft_result[1:(n %/% 2)])^2
  freqs <- (0:(n %/% 2 - 1)) / (n * dt)

  # Find peaks: diel at 1/86400 Hz, annual at 1/(365.25*86400) Hz
  diel_bin <- which.min(abs(freqs - 1/86400))
  annual_bin <- which.min(abs(freqs - 1/(365.25 * 86400)))

  # Diel and annual should be the two dominant peaks
  top_bins <- order(power, decreasing = TRUE)[1:10]
  expect_true(diel_bin %in% top_bins)
  expect_true(annual_bin %in% top_bins)
})


test_that("generate_composite_boundary respects parameters", {

  z_warm <- generate_composite_boundary(T_mu = 20)
  z_cold <- generate_composite_boundary(T_mu = 5)

  expect_gt(mean(zoo::coredata(z_warm)), mean(zoo::coredata(z_cold)))
})


# === S3 constructor input validation =========================================

test_that("thAquifer validates units argument", {
  expect_error(
    thAquifer(0.25, 0.0016, 0.000598, 0.84, 4.186, 3000, 998,
              specificUnits = "not_units"),
    "thUnits"
  )
})


test_that("thHydro validates aquifer argument", {
  expect_error(
    thHydro(hydCond = 10/86400, headGrad = 0.1, aquifer = "not_aquifer"),
    "thAquifer"
  )
})


test_that("thHydro validates unit consistency", {
  aq_default <- thAquifer(0.25, 0.0016, 0.000598, 0.84, 4.186, 3000, 998)
  different_units <- thUnits(L = "ft", t = "hr")

  expect_error(
    thHydro(hydCond = 10/86400, headGrad = 0.1, aquifer = aq_default,
            specificUnits = different_units),
    "Units.*do not match"
  )
})


test_that("thSignal validates unit consistency between hydro and boundary", {
  aq <- thAquifer(0.25, 0.0016, 0.000598, 0.84, 4.186, 3000, 998)
  hy <- thHydro(10/86400, 0.001, 0.1, aq)
  bo_diff <- thBoundary(22, 6, 0, 86400, specificUnits = thUnits(L = "ft"))

  expect_error(thSignal(hy, bo_diff), "must be identical")
})


test_that("thBoundary computes frequency from period", {
  bo <- thBoundary(mean = 22, amplitude = 6, phase = 0, period = 86400)

  expect_equal(bo$frequency, 1/86400)
  expect_true(is.thBoundary(bo))
  expect_true(is.temperheic(bo))
})


test_that("thAquifer computes bulk properties correctly", {
  aq <- thAquifer(porosity = 0.25,
                  thermCond_sed = 0.0016, thermCond_h2o = 0.000598,
                  spHeat_sed = 0.84, spHeat_h2o = 4.186,
                  density_sed = 3000, density_h2o = 998)

  # Bulk thermal conductivity: weighted average
  expected_k <- 0.0016 * 0.75 + 0.000598 * 0.25
  expect_equal(aq$thermCond_bulk, expected_k)

  # Vol heat capacity of water
  expected_vhc_w <- 4.186 * 998
  expect_equal(aq$volHeatCap_h2o, expected_vhc_w)

  # Bulk vol heat capacity: weighted average
  expected_vhc_bulk <- (0.84 * 3000) * 0.75 + expected_vhc_w * 0.25
  expect_equal(aq$volHeatCap_bulk, expected_vhc_bulk)
})

