# stress_test_tier1_spectral.R — Tier 1: Pure spectral and fitting tests
#
# Tests compute_power(), find_peaks(), fit_ols(), fit_fft() with known signals.
# No hydrological complexity — just signal processing correctness.
# Uses both synthetic sinusoids AND real observed signals where spectral
# content is well-characterized (e.g., strong diel in channel 3).
#
# Design: each section is self-contained and can be run independently.
# Summary statistics are collected into a results tibble at the end.
#
# References:
#   Luce et al. (2013) — η parameter, amplitude ratio / phase difference
#   Luce & Tonina (2017) — multi-frequency η consistency
#   Shumway & Stoffer (2017) — spectral analysis fundamentals
#
# Usage: source("scripts/stress_test_tier1_spectral.R") after devtools::load_all()

library(temperheic)
library(zoo)
library(dplyr)
library(purrr)
library(tibble)
Sys.setenv(TZ = "UTC")

# --- Data loading helpers ---
# Robust date parsing: handles both '%Y-%m-%d %H:%M:%S' and '%Y-%m-%d'
parse_dates <- function(x) {
  dt <- as.POSIXct(x, tz = 'UTC', format = '%Y-%m-%d %H:%M:%S')
  missing <- is.na(dt)
  if (any(missing)) {
    dt[missing] <- as.POSIXct(x[missing], tz = 'UTC', format = '%Y-%m-%d')
  }
  dt
}

# Build a clean zoo from CSV: parse dates, remove NAs, deduplicate
load_zoo <- function(filepath, value_col) {
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  date_col <- grep('date|time', names(df), ignore.case = TRUE, value = TRUE)[1]
  dt <- parse_dates(df[[date_col]])
  valid <- !is.na(dt) & !is.na(df[[value_col]])
  idx <- as.numeric(dt[valid])
  z <- zoo(df[[value_col]][valid], order.by = idx)
  if (anyDuplicated(idx)) {
    z <- aggregate(z, by = index(z), FUN = mean, na.rm = TRUE)
  }
  z
}


results <- tibble(
  test_id   = character(),
  test_name = character(),
  metric    = character(),
  expected  = numeric(),
  observed  = numeric(),
  rel_error = numeric(),
  pass      = logical()
)

add_result <- function(id, name, metric, expected, observed, tol = 0.01) {
  re <- if (expected == 0) abs(observed) else abs(observed - expected) / abs(expected)
  results <<- bind_rows(results, tibble(
    test_id = id, test_name = name, metric = metric,
    expected = expected, observed = round(observed, 6),
    rel_error = round(re, 6), pass = re <= tol
  ))
}

# Helper: build an hourly zoo in seconds
make_zoo <- function(vals, n_days, dt = 3600) {
  t_sec <- seq(0, by = dt, length.out = n_days * (86400 / dt))
  if (length(vals) != length(t_sec)) {
    stop("vals length (", length(vals), ") != t_sec length (", length(t_sec), ")")
  }
  zoo(vals, order.by = t_sec)
}

cat("\n========================================\n")
cat("  TIER 1: Pure Spectral & Fitting Tests  \n")
cat("========================================\n\n")


# ============================================================================
# TEST 1.1 — Single cosine recovery
# ============================================================================
cat("--- Test 1.1: Single cosine recovery ---\n")

n_days <- 30
dt <- 3600
t_sec <- seq(0, by = dt, length.out = n_days * 24)
A_true <- 5.0
mean_true <- 15.0
period_true <- 86400

vals_1.1 <- A_true * cos(2 * pi * t_sec / period_true) + mean_true
z_1.1 <- zoo(vals_1.1, order.by = t_sec)

# Spectral analysis
psd_1.1 <- compute_power(z_1.1)
peaks_1.1 <- find_peaks(psd_1.1, prominence = 10, min_separation_hr = 2)

peak_period <- peaks_1.1$period_hr[1]
cat("  Peak period:", round(peak_period, 2), "hr (expected: 24)\n")
add_result("1.1", "Single cosine", "peak_period_hr", 24, peak_period, tol = 0.5/24)

# OLS fit
z_1.1_multi <- zoo(cbind(s1 = vals_1.1), order.by = t_sec)
fit_ols_1.1 <- fit_ols(z_1.1_multi, boundaryMean = mean_true,
                        periodInSeconds = period_true,
                        optimizeRange = c(-0.125, 0.875),
                        nmin = 10, empiricalDataPeriods = n_days)
amp_ols <- attr(fit_ols_1.1, "amplitudes")[1]
cat("  OLS amplitude:", round(amp_ols, 6), "(expected:", A_true, ")\n")
add_result("1.1", "Single cosine", "ols_amplitude", A_true, amp_ols, tol = 1e-10)

# FFT fit
fit_fft_1.1 <- fit_fft(z_1.1_multi, boundaryMean = mean_true,
                        periodInSeconds = period_true,
                        optimizeRange = c(-0.125, 0.875),
                        nmin = 10, empiricalDataPeriods = n_days)
amp_fft <- attr(fit_fft_1.1, "amplitudes")[1]
cat("  FFT amplitude:", round(amp_fft, 6), "(expected:", A_true, ")\n")
add_result("1.1", "Single cosine", "fft_amplitude", A_true, amp_fft, tol = 0.01)

# OLS vs FFT agreement
add_result("1.1", "Single cosine", "ols_fft_amp_agreement",
           amp_ols, amp_fft, tol = 0.01)
cat("\n")


# ============================================================================
# TEST 1.2 — Two-frequency separation (diel + semi-diurnal)
# ============================================================================
cat("--- Test 1.2: Two-frequency separation ---\n")

A_diel <- 5.0
A_semi <- 1.5
vals_1.2 <- A_diel * cos(2 * pi * t_sec / 86400) +
            A_semi * cos(2 * pi * t_sec / 43200) + 15
z_1.2 <- zoo(cbind(s1 = vals_1.2), order.by = t_sec)

psd_1.2 <- compute_power(zoo(vals_1.2, order.by = t_sec))
peaks_1.2 <- find_peaks(psd_1.2, prominence = 10, min_separation_hr = 2)

cat("  Peaks found:", nrow(peaks_1.2), "(expected: 2)\n")
top2 <- peaks_1.2[1:min(2, nrow(peaks_1.2)), ]
periods_sorted <- sort(top2$period_hr)
cat("  Peak periods:", round(periods_sorted, 1), "hr (expected: 12, 24)\n")
add_result("1.2", "Two-frequency", "peak1_period_hr", 12, periods_sorted[1], tol = 0.5/12)
add_result("1.2", "Two-frequency", "peak2_period_hr", 24, periods_sorted[2], tol = 0.5/24)

# OLS at diel — should recover 5.0 without contamination from 12h
fit_diel <- fit_ols(z_1.2, boundaryMean = 15, periodInSeconds = 86400,
                    optimizeRange = c(-0.125, 0.875), nmin = 10,
                    empiricalDataPeriods = n_days)
amp_diel_ols <- attr(fit_diel, "amplitudes")[1]
cat("  OLS diel amplitude:", round(amp_diel_ols, 4), "(expected:", A_diel, ")\n")
add_result("1.2", "Two-frequency", "ols_diel_amp", A_diel, amp_diel_ols, tol = 0.01)

# OLS at semi-diurnal — should recover 1.5
fit_semi <- fit_ols(z_1.2, boundaryMean = 15, periodInSeconds = 43200,
                    optimizeRange = c(-0.125, 0.875), nmin = 10,
                    empiricalDataPeriods = n_days * 2)
amp_semi_ols <- attr(fit_semi, "amplitudes")[1]
cat("  OLS semi-diurnal amplitude:", round(amp_semi_ols, 4), "(expected:", A_semi, ")\n")
add_result("1.2", "Two-frequency", "ols_semi_amp", A_semi, amp_semi_ols, tol = 0.01)
cat("\n")


# ============================================================================
# TEST 1.3 — Diel + annual (multi-scale)
# ============================================================================
cat("--- Test 1.3: Diel + annual (2 years) ---\n")

n_days_long <- 730
t_long <- seq(0, by = 3600, length.out = n_days_long * 24)
A_diel_1.3 <- 3.0
A_ann_1.3  <- 10.0
vals_1.3 <- A_diel_1.3 * cos(2 * pi * t_long / 86400) +
            A_ann_1.3  * cos(2 * pi * t_long / (365.25 * 86400)) + 12
z_1.3 <- zoo(vals_1.3, order.by = t_long)

# Spectrum should show both peaks
psd_1.3 <- compute_power(z_1.3)
peaks_1.3 <- find_peaks(psd_1.3, prominence = 20, min_separation_hr = 10)
cat("  Peaks found:", nrow(peaks_1.3), "\n")
if (nrow(peaks_1.3) >= 2) {
  cat("  Top 2 periods:", round(peaks_1.3$period_hr[1:2], 1), "hr\n")
  cat("  Top 2 periods:", round(peaks_1.3$period_day[1:2], 1), "days\n")
}

# Detrend and recover diel
z_1.3_detrended <- detrend_temperature(z_1.3, method = "loess", span = 0.02)
z_1.3_det_multi <- zoo(cbind(s1 = coredata(z_1.3_detrended)), order.by = index(z_1.3))
fit_1.3_diel <- fit_ols(z_1.3_det_multi, boundaryMean = 0,
                         periodInSeconds = 86400,
                         optimizeRange = c(-0.125, 0.875),
                         nmin = 10, empiricalDataPeriods = n_days_long)
amp_diel_detr <- attr(fit_1.3_diel, "amplitudes")[1]
cat("  Detrended diel amplitude:", round(amp_diel_detr, 4), "(expected:", A_diel_1.3, ")\n")
add_result("1.3", "Diel+annual", "detrended_diel_amp", A_diel_1.3, amp_diel_detr, tol = 0.05)

# Annual fit on full record
z_1.3_multi <- zoo(cbind(s1 = vals_1.3), order.by = t_long)
fit_1.3_ann <- fit_ols(z_1.3_multi, boundaryMean = 12,
                        periodInSeconds = 365.25 * 86400,
                        optimizeRange = c(-0.125, 0.875),
                        nmin = 10, empiricalDataPeriods = 2)
amp_ann_full <- attr(fit_1.3_ann, "amplitudes")[1]
cat("  Full-record annual amplitude:", round(amp_ann_full, 4), "(expected:", A_ann_1.3, ")\n")
add_result("1.3", "Diel+annual", "annual_amp", A_ann_1.3, amp_ann_full, tol = 0.01)

# Windowed diel stability (7-day windows on detrended signal)
windows_1.3 <- window_temperature(z_1.3_detrended, width = 7 * 86400, step = 86400)
windowed_amps <- map_dbl(windows_1.3$data, function(w) {
  wz <- zoo(cbind(s1 = coredata(w)), order.by = index(w))
  f <- fit_ols(wz, boundaryMean = 0, periodInSeconds = 86400,
               optimizeRange = c(-0.125, 0.875), nmin = 10,
               empiricalDataPeriods = 7)
  attr(f, "amplitudes")[1]
})
cat("  Windowed diel amplitude: median =", round(median(windowed_amps), 4),
    ", CV =", round(sd(windowed_amps)/mean(windowed_amps), 4), "\n")
add_result("1.3", "Diel+annual", "windowed_diel_cv", 0, 
           sd(windowed_amps)/mean(windowed_amps), tol = 0.05)
cat("\n")


# ============================================================================
# TEST 1.4 — Non-sinusoidal diel (harmonics)
# ============================================================================
cat("--- Test 1.4: Non-sinusoidal diel (3 harmonics) ---\n")

A1 <- 4.0; A2 <- 1.5; A3 <- 0.6
vals_1.4 <- A1 * cos(2 * pi * t_sec / 86400) +
            A2 * cos(2 * pi * t_sec / 43200) +
            A3 * cos(2 * pi * t_sec / 28800) + 15
z_1.4 <- zoo(vals_1.4, order.by = t_sec)

psd_1.4 <- compute_power(z_1.4)
peaks_1.4 <- find_peaks(psd_1.4, prominence = 5, min_separation_hr = 2)
cat("  Peaks found:", nrow(peaks_1.4), "(expected: 3)\n")
add_result("1.4", "Non-sinusoidal", "n_peaks", 3, nrow(peaks_1.4), tol = 0)

if (nrow(peaks_1.4) >= 3) {
  periods_1.4 <- sort(peaks_1.4$period_hr[1:3])
  cat("  Periods:", round(periods_1.4, 1), "hr (expected: 8, 12, 24)\n")
  add_result("1.4", "Non-sinusoidal", "harmonic_8h", 8, periods_1.4[1], tol = 0.5/8)
  add_result("1.4", "Non-sinusoidal", "harmonic_12h", 12, periods_1.4[2], tol = 0.5/12)
  add_result("1.4", "Non-sinusoidal", "harmonic_24h", 24, periods_1.4[3], tol = 0.5/24)
}

# OLS at each harmonic
z_1.4_multi <- zoo(cbind(s1 = vals_1.4), order.by = t_sec)
for (info in list(c(86400, A1, "24h"), c(43200, A2, "12h"), c(28800, A3, "8h"))) {
  per <- as.numeric(info[1]); amp_true <- as.numeric(info[2]); label <- info[3]
  f <- fit_ols(z_1.4_multi, boundaryMean = 15, periodInSeconds = per,
               optimizeRange = c(-0.125, 0.875), nmin = 10,
               empiricalDataPeriods = n_days * 86400 / per)
  amp_recovered <- attr(f, "amplitudes")[1]
  cat("  OLS", label, "amplitude:", round(amp_recovered, 4), "(expected:", amp_true, ")\n")
  add_result("1.4", "Non-sinusoidal", paste0("ols_amp_", label), amp_true, 
             amp_recovered, tol = 0.01)
}
cat("\n")


# ============================================================================
# TEST 1.5 — Record length sensitivity
# ============================================================================
cat("--- Test 1.5: Record length sensitivity ---\n")

for (nd in c(3, 7, 14, 30, 90)) {
  t_tmp <- seq(0, by = 3600, length.out = nd * 24)
  v_tmp <- 5 * cos(2 * pi * t_tmp / 86400) + 
           1.5 * cos(2 * pi * t_tmp / 43200) + 15
  z_tmp <- zoo(cbind(s1 = v_tmp), order.by = t_tmp)
  
  f_tmp <- fit_ols(z_tmp, boundaryMean = 15, periodInSeconds = 86400,
                   optimizeRange = c(-0.125, 0.875), nmin = 10,
                   empiricalDataPeriods = nd)
  amp_tmp <- attr(f_tmp, "amplitudes")[1]
  err_tmp <- abs(amp_tmp - 5.0) / 5.0
  cat("  ", nd, "days: OLS diel amp =", round(amp_tmp, 4), 
      ", rel error =", round(err_tmp, 6), "\n")
  add_result("1.5", "Record length", paste0("ols_amp_", nd, "d"),
             5.0, amp_tmp, tol = 0.01)
}
cat("\n")


# ============================================================================
# TEST 1.6 — Phase recovery (two-sensor pairwise)
# ============================================================================
cat("--- Test 1.6: Phase recovery (two sensors) ---\n")

A1_true <- 5.0; A2_true <- 3.0
phi1_true <- 0.0; phi2_true <- 0.5  # radians
ar_true <- A2_true / A1_true
dphi_true <- phi2_true - phi1_true
eta_true <- -log(ar_true) / dphi_true

s1 <- A1_true * cos(2 * pi * t_sec / 86400 - phi1_true) + 15
s2 <- A2_true * cos(2 * pi * t_sec / 86400 - phi2_true) + 15
z_1.6 <- zoo(cbind(sensor1 = s1, sensor2 = s2), order.by = t_sec)

fit_1.6 <- fit_ols(z_1.6, boundaryMean = 15, periodInSeconds = 86400,
                    optimizeRange = c(-0.125, 0.875), nmin = 10,
                    empiricalDataPeriods = n_days)

amps_1.6 <- attr(fit_1.6, "amplitudes")
phases_1.6 <- attr(fit_1.6, "phases")
ar_recovered <- amps_1.6[2] / amps_1.6[1]
# Phase difference: the fit_ols pairwise output
# deltaPhaseRadians is a length-4 vector for 2 sensors (2x2 matrix flattened)
dphi_recovered <- fit_1.6$deltaPhaseRadians
# The (1,2) element: sensor2 - sensor1
dphi_12 <- dphi_recovered[3]  # expand.grid ordering: (1,1), (2,1), (1,2), (2,2)

cat("  Amp ratio:", round(ar_recovered, 6), "(expected:", ar_true, ")\n")
cat("  Phase diff:", round(dphi_12, 6), "rad (expected:", dphi_true, ")\n")
add_result("1.6", "Phase recovery", "amp_ratio", ar_true, ar_recovered, tol = 1e-10)
add_result("1.6", "Phase recovery", "phase_diff", dphi_true, dphi_12, tol = 1e-6)

# Compute η from recovered values
eta_recovered <- -log(ar_recovered) / dphi_12
cat("  η:", round(eta_recovered, 6), "(expected:", round(eta_true, 6), ")\n")
add_result("1.6", "Phase recovery", "eta", eta_true, eta_recovered, tol = 1e-6)
cat("\n")


# ============================================================================
# TEST 1.7 — Real data: channel 3 spectral characterization
# ============================================================================
cat("--- Test 1.7: Channel 3 observed spectrum ---\n")

data_dir <- system.file("extdata", package = "temperheic")
ch3_zoo <- load_zoo(file.path(data_dir, "meacham_channel3_hourly.csv"), "channel_3_temp")

psd_ch3 <- compute_power(ch3_zoo)
peaks_ch3 <- find_peaks(psd_ch3, prominence = 20, min_separation_hr = 3)

# Should find diel peak
diel_peak <- peaks_ch3[which.min(abs(peaks_ch3$period_hr - 24)), ]
has_diel <- nrow(diel_peak) > 0 && abs(diel_peak$period_hr - 24) < 1
cat("  Diel peak found:", has_diel, "\n")
if (has_diel) cat("  Diel peak period:", round(diel_peak$period_hr, 2), "hr\n")
add_result("1.7", "Channel 3 spectrum", "diel_detected", 1, as.numeric(has_diel), tol = 0)

# Should find semi-diurnal
semi_peak <- peaks_ch3[which.min(abs(peaks_ch3$period_hr - 12)), ]
has_semi <- nrow(semi_peak) > 0 && abs(semi_peak$period_hr - 12) < 1
cat("  Semi-diurnal peak found:", has_semi, "\n")
add_result("1.7", "Channel 3 spectrum", "semi_diurnal_detected", 1, 
           as.numeric(has_semi), tol = 0)

# OLS diel amplitude should be in 3-8°C range (physically reasonable)
ch3_multi <- zoo(cbind(ch3 = coredata(ch3_zoo)), order.by = index(ch3_zoo))
fit_ch3 <- fit_ols(ch3_multi, boundaryMean = mean(coredata(ch3_zoo)),
                   periodInSeconds = 86400,
                   optimizeRange = c(-0.125, 0.875), nmin = 10,
                   empiricalDataPeriods = round(diff(range(index(ch3_zoo))) / 86400))
amp_ch3 <- attr(fit_ch3, "amplitudes")[1]
cat("  Full-record OLS diel amplitude:", round(amp_ch3, 2), "°C\n")
reasonable_amp <- amp_ch3 >= 1.0 & amp_ch3 <= 10.0
add_result("1.7", "Channel 3 spectrum", "diel_amp_reasonable", 1, 
           as.numeric(reasonable_amp), tol = 0)
cat("\n")


# ============================================================================
# TEST 1.8 — Real data: spring 6 null case
# ============================================================================
cat("--- Test 1.8: Spring 6 null case (no diel expected) ---\n")

sp6_zoo <- load_zoo(file.path(data_dir, "meacham_spring6_hourly.csv"), "spring_6_temp")

psd_sp6 <- compute_power(sp6_zoo)
peaks_sp6 <- find_peaks(psd_sp6, prominence = 50, min_separation_hr = 3)

# Should NOT find a prominent diel peak
diel_cand <- peaks_sp6[abs(peaks_sp6$period_hr - 24) < 1, ]
has_diel_sp6 <- nrow(diel_cand) > 0
cat("  Diel peak detected (prominence>50):", has_diel_sp6, "(expected: FALSE)\n")
add_result("1.8", "Spring 6 null", "no_diel_detected", 0, 
           as.numeric(has_diel_sp6), tol = 0)

# OLS diel amplitude should be very small (<0.5°C)
sp6_multi <- zoo(cbind(sp6 = coredata(sp6_zoo)), order.by = index(sp6_zoo))
fit_sp6 <- fit_ols(sp6_multi, boundaryMean = mean(coredata(sp6_zoo)),
                   periodInSeconds = 86400,
                   optimizeRange = c(-0.125, 0.875), nmin = 10,
                   empiricalDataPeriods = round(diff(range(index(sp6_zoo))) / 86400))
amp_sp6 <- attr(fit_sp6, "amplitudes")[1]
cat("  Spring 6 diel amplitude:", round(amp_sp6, 4), "°C (should be < 0.5)\n")
small_amp <- amp_sp6 < 0.5
add_result("1.8", "Spring 6 null", "diel_amp_small", 1, as.numeric(small_amp), tol = 0)
cat("\n")


# ============================================================================
# TEST 1.9 — Real data: multi-site spectral comparison
# ============================================================================
cat("--- Test 1.9: Multi-site spectral comparison ---\n")
cat("  Comparing diel PSD across signal types\n")

# Air temperature
air_zoo <- load_zoo(file.path(data_dir, "UP211_air_temp_hourly_2012_2014.csv"), "air_temp")

# Well 20
w20_zoo <- load_zoo(file.path(data_dir, "meacham_well20_hourly.csv"), "well_20_temp")

# Well 10
w10_zoo <- load_zoo(file.path(data_dir, "meacham_well10_hourly.csv"), "well_10_temp")

# Compute diel PSD for each
get_diel_psd <- function(z, label) {
  psd <- compute_power(z)
  diel_idx <- which.min(abs(psd$period_hr - 24))
  tibble(label = label, diel_psd = psd$psd[diel_idx], period_hr = psd$period_hr[diel_idx])
}

diel_comparison <- bind_rows(
  get_diel_psd(air_zoo, "Air"),
  get_diel_psd(ch3_zoo, "Channel 3"),
  get_diel_psd(w20_zoo, "Well 20"),
  get_diel_psd(w10_zoo, "Well 10"),
  get_diel_psd(sp6_zoo, "Spring 6")
)

cat("  Diel PSD by site:\n")
for (i in seq_len(nrow(diel_comparison))) {
  cat("    ", diel_comparison$label[i], ":", round(diel_comparison$diel_psd[i], 1), "\n")
}

# Physical ordering: Air > Stream > Well 20 > Well 10 ≈ Spring 6
air_psd <- diel_comparison$diel_psd[diel_comparison$label == "Air"]
ch3_psd <- diel_comparison$diel_psd[diel_comparison$label == "Channel 3"]
w20_psd <- diel_comparison$diel_psd[diel_comparison$label == "Well 20"]
w10_psd <- diel_comparison$diel_psd[diel_comparison$label == "Well 10"]
sp6_psd <- diel_comparison$diel_psd[diel_comparison$label == "Spring 6"]

ordering_correct <- (air_psd > ch3_psd) & (ch3_psd > w20_psd) & (w20_psd > w10_psd)
cat("  Physical ordering (Air > Stream > Well20 > Well10):", ordering_correct, "\n")
add_result("1.9", "Multi-site comparison", "diel_ordering", 1, 
           as.numeric(ordering_correct), tol = 0)
cat("\n")


# ============================================================================
# TEST 1.10 — Windowed diel amplitude seasonal envelope (real data)
# ============================================================================
cat("--- Test 1.10: Windowed diel amplitude seasonal envelope ---\n")

# Window channel 3 into 7-day windows, extract diel amplitude per window
# Should show seasonal modulation: larger in summer, smaller in winter
ch3_windows <- window_temperature(ch3_zoo, width = 7 * 86400, step = 86400)
cat("  Windows:", nrow(ch3_windows), "\n")

ch3_window_amps <- map_dbl(ch3_windows$data, function(w) {
  # Detrend within window
  w_det <- detrend_temperature(w, method = "polynomial", degree = 1)
  wz <- zoo(cbind(s1 = coredata(w_det)), order.by = index(w))
  tryCatch({
    f <- fit_ols(wz, boundaryMean = 0, periodInSeconds = 86400,
                 optimizeRange = c(-0.125, 0.875), nmin = 10,
                 empiricalDataPeriods = 7)
    attr(f, "amplitudes")[1]
  }, error = function(e) NA_real_)
})

valid_amps <- ch3_window_amps[!is.na(ch3_window_amps)]
cat("  Valid windows:", length(valid_amps), "/", nrow(ch3_windows), "\n")
cat("  Amplitude range:", round(range(valid_amps), 2), "°C\n")
cat("  Amplitude median:", round(median(valid_amps), 2), "°C\n")

# Summer amps should be larger than winter amps
# Window centers as dates for seasonal check
window_dates <- as.POSIXct(ch3_windows$window_center, origin = "1970-01-01", tz = "UTC")
months <- as.numeric(format(window_dates, "%m"))
summer <- months %in% c(6, 7, 8)
winter <- months %in% c(12, 1, 2)

summer_median <- median(ch3_window_amps[summer], na.rm = TRUE)
winter_median <- median(ch3_window_amps[winter], na.rm = TRUE)
cat("  Summer median:", round(summer_median, 2), "°C\n")
cat("  Winter median:", round(winter_median, 2), "°C\n")
seasonal_modulation <- summer_median > winter_median
cat("  Summer > Winter:", seasonal_modulation, "\n")
add_result("1.10", "Seasonal envelope", "summer_gt_winter", 1,
           as.numeric(seasonal_modulation), tol = 0)
cat("\n")


# ============================================================================
# SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("  TIER 1 RESULTS SUMMARY\n")
cat("========================================\n\n")

cat("Total tests:", nrow(results), "\n")
cat("Passed:", sum(results$pass), "\n")
cat("Failed:", sum(!results$pass), "\n\n")

if (any(!results$pass)) {
  cat("FAILURES:\n")
  print(results[!results$pass, ], n = Inf)
} else {
  cat("All tests passed.\n")
}

cat("\nFull results:\n")
print(results, n = Inf)
