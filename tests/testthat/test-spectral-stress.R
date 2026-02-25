# tests/testthat/test-spectral-stress.R
#
# Stress tests for spectral analysis pipeline (Tier 1).
# Complements test-compute_power.R (structure/basic correctness) and
# test-round-trip.R (forward-inverse eta recovery) with harder scenarios:
#   - Multi-scale signal separation (diel + annual)
#   - Harmonic amplitude recovery via find_peaks()
#   - find_peaks() output structure and input validation
#   - Windowed amplitude recovery (Hann correction factor)

# --- Shared helper ---
make_hourly_zoo <- function(n_days, components, mean_temp = 10) {
  # components: list of list(amplitude, period_s, phase = 0)
  t_sec <- seq(0, by = 3600, length.out = n_days * 24)
  vals <- rep(mean_temp, length(t_sec))
  for (comp in components) {
    ph <- if (!is.null(comp$phase)) comp$phase else 0
    vals <- vals + comp$amplitude * cos(2 * pi * t_sec / comp$period_s - ph)
  }
  zoo::zoo(vals, order.by = t_sec)
}


# === find_peaks: output structure =============================================

test_that("find_peaks returns expected columns", {
  z <- make_hourly_zoo(30, list(list(amplitude = 5, period_s = 86400)))
  psd <- compute_power(z)
  peaks <- find_peaks(psd, n_peaks = 3)

  expect_s3_class(peaks, "tbl_df")
  expect_true(all(c("frequency_hz", "period_s", "period_hr",
                     "period_day", "psd", "amplitude") %in% names(peaks)))
})

test_that("find_peaks returns peaks sorted by PSD descending", {
  z <- make_hourly_zoo(30, list(
    list(amplitude = 5, period_s = 86400),
    list(amplitude = 2, period_s = 43200)
  ))
  psd <- compute_power(z)
  peaks <- find_peaks(psd, n_peaks = 5)

  expect_equal(peaks$psd, sort(peaks$psd, decreasing = TRUE))
})

test_that("find_peaks n_peaks limits output rows", {
  z <- make_hourly_zoo(30, list(
    list(amplitude = 5, period_s = 86400),
    list(amplitude = 2, period_s = 43200),
    list(amplitude = 1, period_s = 28800)
  ))
  psd <- compute_power(z)
  peaks_2 <- find_peaks(psd, n_peaks = 2)
  expect_lte(nrow(peaks_2), 2)
})


# === find_peaks: amplitude recovery ==========================================

test_that("find_peaks recovers exact amplitude for single cosine", {
  z <- make_hourly_zoo(30, list(list(amplitude = 5, period_s = 86400)))
  psd <- compute_power(z)
  peaks <- find_peaks(psd, n_peaks = 1)

  expect_equal(peaks$amplitude[1], 5, tolerance = 0.05)
})

test_that("find_peaks recovers amplitudes for three harmonics", {
  # Non-sinusoidal diel: fundamental + 2nd + 3rd harmonic
  # (Luce 2017 multi-frequency basis)
  z <- make_hourly_zoo(30, list(
    list(amplitude = 5.0, period_s = 86400),
    list(amplitude = 1.5, period_s = 43200),
    list(amplitude = 0.8, period_s = 28800)
  ))
  psd <- compute_power(z)
  peaks <- find_peaks(psd, n_peaks = 5)

  expect_gte(nrow(peaks), 3)

  p24 <- peaks[which.min(abs(peaks$period_hr - 24)), ]
  p12 <- peaks[which.min(abs(peaks$period_hr - 12)), ]
  p8  <- peaks[which.min(abs(peaks$period_hr - 8)),  ]

  expect_equal(p24$amplitude, 5.0, tolerance = 0.1)
  expect_equal(p12$amplitude, 1.5, tolerance = 0.1)
  expect_equal(p8$amplitude,  0.8, tolerance = 0.1)
})


# === Multi-scale: diel + annual ===============================================

test_that("compute_power separates diel and annual signals over 2-year record", {
  z <- make_hourly_zoo(730, list(
    list(amplitude = 8, period_s = 365.25 * 86400),
    list(amplitude = 3, period_s = 86400)
  ))
  psd <- compute_power(z)
  peaks <- find_peaks(psd, n_peaks = 5)

  expect_gte(nrow(peaks), 2)
  annual_peak <- peaks[which.min(abs(peaks$period_day - 365)), ]
  diel_peak   <- peaks[which.min(abs(peaks$period_hr - 24)),   ]

  expect_equal(annual_peak$period_day, 365, tolerance = 5)
  expect_equal(diel_peak$period_hr, 24, tolerance = 0.5)

  expect_gt(annual_peak$psd, diel_peak$psd)

  expect_equal(annual_peak$amplitude, 8, tolerance = 0.5)
  expect_equal(diel_peak$amplitude,   3, tolerance = 0.2)
})


# === Windowing: Hann correction ===============================================

test_that("Hann window preserves peak location but reduces PSD", {
  z <- make_hourly_zoo(30, list(list(amplitude = 5, period_s = 86400)))
  psd_rect <- compute_power(z, window = "none")
  psd_hann <- compute_power(z, window = "hann")

  peak_rect <- psd_rect[which.max(psd_rect$psd), ]
  peak_hann <- psd_hann[which.max(psd_hann$psd), ]

  expect_equal(peak_rect$period_hr, peak_hann$period_hr)
  expect_lt(peak_hann$psd, peak_rect$psd)
})


# === Input validation =========================================================

test_that("find_peaks rejects non-tibble input", {
  expect_error(find_peaks(data.frame(x = 1:10)))
})

test_that("find_peaks handles white noise gracefully", {
  set.seed(42)
  t_sec <- seq(0, by = 3600, length.out = 720)
  z <- zoo::zoo(rnorm(720), order.by = t_sec)
  psd <- compute_power(z)
  peaks <- find_peaks(psd, n_peaks = 3)

  expect_lte(nrow(peaks), 3)
  expect_true(all(peaks$psd > 0))
})
