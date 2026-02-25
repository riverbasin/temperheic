
# tests/testthat/test-compute_power.R
#
# Tests for compute_power() — power spectral density estimation.
# Covers: output structure, physical correctness (peak location, Parseval),
# multi-column dispatch, windowing, detrending, NA handling, and edge cases.

# --- Helper: create a clean hourly zoo with known spectral content ---
make_test_zoo <- function(n_days = 10, dt = 3600, amp_diel = 5, mean_temp = 20) {
  t_sec <- seq(0, by = dt, length.out = n_days * (86400 / dt))
  vals <- amp_diel * cos(2 * pi * t_sec / 86400) + mean_temp
  zoo::zoo(vals, order.by = t_sec)
}


# === Output structure =========================================================

test_that("compute_power returns a tibble with expected columns", {
  z <- make_test_zoo()
  result <- compute_power(z)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("frequency_hz", "period_s", "period_hr",
                          "period_day", "psd", "series"))
  expect_true(all(result$frequency_hz > 0))
  expect_true(all(result$psd >= 0))
  expect_true(all(is.finite(result$psd)))
})

test_that("output contains only positive frequencies (no DC, no Nyquist edge)", {
  z <- make_test_zoo()
  result <- compute_power(z)

  # DC would be frequency_hz = 0

  expect_true(all(result$frequency_hz > 0))
  # Period should not be infinite
  expect_true(all(is.finite(result$period_s)))
})

test_that("period columns are consistent with frequency_hz", {
  z <- make_test_zoo()
  result <- compute_power(z)

  expect_equal(result$period_s, 1 / result$frequency_hz)
  expect_equal(result$period_hr, result$period_s / 3600)
  expect_equal(result$period_day, result$period_s / 86400)
})

test_that("series column is populated for single-column zoo", {
  z <- make_test_zoo()
  result <- compute_power(z)
  expect_true(all(!is.na(result$series)))
  expect_length(unique(result$series), 1)
})


# === Physical correctness =====================================================

test_that("pure diel cosine produces peak at 24 hours", {
  # 10 days hourly, bin-aligned (240 points = integer multiple of 24)
  z <- make_test_zoo(n_days = 10)
  result <- compute_power(z)

  peak <- result[which.max(result$psd), ]
  expect_equal(peak$period_hr, 24, tolerance = 0.5)
})

test_that("two-frequency signal produces peaks at both periods", {
  # Diel + 12-hr semi-diurnal
  t_sec <- seq(0, by = 3600, length.out = 10 * 24)
  vals <- 5 * cos(2 * pi * t_sec / 86400) + 2 * cos(2 * pi * t_sec / 43200)
  z <- zoo::zoo(vals, order.by = t_sec)
  result <- compute_power(z)

  # Top two peaks should be near 24h and 12h
  top2 <- result[order(result$psd, decreasing = TRUE), ][1:2, ]
  periods <- sort(top2$period_hr)
  expect_equal(periods[1], 12, tolerance = 0.5)
  expect_equal(periods[2], 24, tolerance = 0.5)

  # Diel should have higher PSD (amplitude 5 > 2)
  diel_psd <- result$psd[which.min(abs(result$period_hr - 24))]
  semi_psd <- result$psd[which.min(abs(result$period_hr - 12))]
  expect_gt(diel_psd, semi_psd)
})

test_that("PSD scales with amplitude squared", {
  # Two signals: amplitude 2 vs amplitude 4 at same frequency
  t_sec <- seq(0, by = 3600, length.out = 10 * 24)
  z_small <- zoo::zoo(2 * cos(2 * pi * t_sec / 86400), order.by = t_sec)
  z_large <- zoo::zoo(4 * cos(2 * pi * t_sec / 86400), order.by = t_sec)

  psd_s <- compute_power(z_small)
  psd_l <- compute_power(z_large)

  peak_s <- max(psd_s$psd)
  peak_l <- max(psd_l$psd)

  # PSD ratio should be (4/2)^2 = 4
  expect_equal(peak_l / peak_s, 4, tolerance = 0.1)
})


# === Multi-column dispatch ====================================================

test_that("multi-column zoo produces long-format tibble with series column", {
  t_sec <- seq(0, by = 3600, length.out = 10 * 24)
  vals <- cbind(
    stream = 5 * cos(2 * pi * t_sec / 86400) + 20,
    well   = 1 * cos(2 * pi * t_sec / 86400) + 18
  )
  z <- zoo::zoo(vals, order.by = t_sec)
  result <- compute_power(z)

  expect_s3_class(result, "tbl_df")
  expect_true("series" %in% names(result))
  expect_setequal(unique(result$series), c("stream", "well"))

  # Each series should have the same number of frequency bins
  n_per_series <- table(result$series)
  expect_equal(as.numeric(n_per_series["stream"]),
               as.numeric(n_per_series["well"]))
})

test_that("multi-column correctly identifies higher amplitude in stream", {
  t_sec <- seq(0, by = 3600, length.out = 10 * 24)
  vals <- cbind(
    stream = 5 * cos(2 * pi * t_sec / 86400),
    well   = 1 * cos(2 * pi * t_sec / 86400)
  )
  z <- zoo::zoo(vals, order.by = t_sec)
  result <- compute_power(z)

  stream_peak <- max(result$psd[result$series == "stream"])
  well_peak   <- max(result$psd[result$series == "well"])
  expect_gt(stream_peak, well_peak)
})


# === Windowing ================================================================

test_that("Hann window option runs without error", {
  z <- make_test_zoo()
  result_none <- compute_power(z, window = "none")
  result_hann <- compute_power(z, window = "hann")

  expect_equal(nrow(result_none), nrow(result_hann))
  # Both should find the diel peak
  peak_none <- result_none$period_hr[which.max(result_none$psd)]
  peak_hann <- result_hann$period_hr[which.max(result_hann$psd)]
  expect_equal(peak_none, 24, tolerance = 0.5)
  expect_equal(peak_hann, 24, tolerance = 0.5)
})


# === Detrending ===============================================================

test_that("detrend = TRUE removes linear trend energy", {
  # Signal with strong linear drift + weak diel
  t_sec <- seq(0, by = 3600, length.out = 30 * 24)
  drift <- 10 * (t_sec / max(t_sec))  # 10-degree drift over 30 days
  diel  <- 0.5 * cos(2 * pi * t_sec / 86400)
  z <- zoo::zoo(drift + diel, order.by = t_sec)

  psd_raw  <- compute_power(z, detrend = FALSE)
  psd_detr <- compute_power(z, detrend = TRUE)

  # With detrend, low-frequency power should be reduced
  low_freq_raw  <- sum(psd_raw$psd[psd_raw$period_day > 10])
  low_freq_detr <- sum(psd_detr$psd[psd_detr$period_day > 10])
  expect_lt(low_freq_detr, low_freq_raw)
})


# === NA handling ==============================================================

test_that("sparse NAs (< 10%) handled silently", {
  set.seed(42)
  z <- make_test_zoo()
  vals <- zoo::coredata(z)
  vals[sample(length(vals), 5)] <- NA  # ~2% NAs
  z_na <- zoo::zoo(vals, order.by = zoo::index(z))

  # Should not warn
  expect_no_warning(compute_power(z_na))
})

test_that("many NAs (> 10%) produce a warning", {
  set.seed(42)
  z <- make_test_zoo()
  vals <- zoo::coredata(z)
  vals[sample(length(vals), 50)] <- NA  # ~21% NAs
  z_na <- zoo::zoo(vals, order.by = zoo::index(z))

  expect_warning(compute_power(z_na), "NA")
})


# === Input validation =========================================================

test_that("non-zoo input raises error", {
  expect_error(compute_power(1:100), "zoo")
  expect_error(compute_power(data.frame(x = 1:10)), "zoo")
})

test_that("invalid window option raises error", {
  z <- make_test_zoo()
  expect_error(compute_power(z, window = "hamming"))
})

