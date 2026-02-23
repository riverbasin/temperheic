
# testingObservedFits.R — Fitting functions applied to real temperature time series
#
# Purpose:
#   Exercise temperheic's fitting and signal processing tools on observed
#   Meacham Creek data. Companion to testingThSeries.R (which uses synthetics).
#   Covers the full spectrum of signal types:
#     - Stream temperature: strong diel + seasonal (the gold standard)
#     - Air temperature: large diel amplitude, noisy (weather), seasonal
#     - Shallow well (well 20): detectable diel signal, attenuated
#     - Deep wells (well 21, well 10): seasonal only, diel fully damped
#     - Spring (offchannel spring 6): seasonal only, no diel, minimal range
#   Ends with the full inverse solution (thObservedSeries) on the 3-site
#   annual analysis for flux estimation.
#
# Data sources (all in inst/extdata/):
#   meacham_3site_hourly.csv  — channel 3 (0m), well 21 (175m), well 10 (955m)
#   UP211_air_temp_hourly_2012_2014.csv — SNOTEL air temp, hourly
#   meacham_well20_hourly.csv — shallow floodplain well with diel signal
#   meacham_spring6_hourly.csv — offchannel spring (seasonal only)
#
# References:
#   Luce et al. (2013) — analytical solutions, eta parameter, inverse model
#   Vogt et al. (2012) — x_d / v_d formulation of the ADE solution
#
# Usage: source interactively in RStudio, or run sections as needed.

library(temperheic)
library(zoo)
library(dplyr)
library(purrr)
Sys.setenv(TZ = "UTC")


# === 1. Load data =============================================================
#
# All CSVs are lightweight extracts from the meachamCreekRestoration.sqlite
# database. Timestamps are UTC. See the extraction code in the session notes
# for provenance.

data_dir <- system.file("extdata", package = "temperheic")

# 3-site data for annual inverse solution
three_site <- read.csv(file.path(data_dir, "meacham_3site_hourly.csv")) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))

# Air temperature (SNOTEL station UP211, Bonifer)
air_raw <- read.csv(file.path(data_dir, "UP211_air_temp_hourly_2012_2014.csv")) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))

# Well 20 — shallow floodplain well with detectable diel signal
well20_raw <- read.csv(file.path(data_dir, "meacham_well20_hourly.csv")) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))

# Offchannel spring 6 — seasonal only
spring6_raw <- read.csv(file.path(data_dir, "meacham_spring6_hourly.csv")) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))

# Full individual records (longer than the 3-site concurrent overlap)
ch3_full <- read.csv(file.path(data_dir, "meacham_channel3_hourly.csv")) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))

w10_full <- read.csv(file.path(data_dir, "meacham_well10_hourly.csv")) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")) %>%
  mutate(well_10_temp = ifelse(well_10_temp == 0, NA, well_10_temp))  # instrument zeros -> NA

sp6_full <- read.csv(file.path(data_dir, "meacham_spring6_hourly.csv")) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))

# Build zoo objects (seconds index for consistency with temperheic convention)
to_zoo_seconds <- function(df, value_col) {
  idx <- as.numeric(df$date_time)
  z <- suppressWarnings(zoo(df[[value_col]], order.by = idx))
  # Deduplicate: average values at identical timestamps
  if (anyDuplicated(idx)) {
    z <- aggregate(z, by = index(z), FUN = mean, na.rm = TRUE)
  }
  z
}

stream_zoo  <- to_zoo_seconds(three_site, "channel_3_0m")
well21_zoo  <- to_zoo_seconds(three_site, "well_21_175m")
well10_zoo  <- to_zoo_seconds(three_site, "well_10_955m")
air_zoo     <- to_zoo_seconds(air_raw, "air_temp")
well20_zoo  <- to_zoo_seconds(well20_raw, "well_20_temp")
spring6_zoo <- to_zoo_seconds(spring6_raw, "spring_6_temp")

# Full individual records as zoo (use for spectral analysis & single-site fits)
stream_full_zoo  <- to_zoo_seconds(ch3_full, "channel_3_temp")
well10_full_zoo  <- to_zoo_seconds(w10_full, "well_10_temp")
spring6_full_zoo <- to_zoo_seconds(sp6_full, "spring_6_temp")

cat("=== Data loaded ===\n")
cat("  Stream (channel 3):", length(stream_zoo), "obs,",
    round(diff(range(index(stream_zoo))) / 86400), "days\n")
cat("  Air (UP211):       ", length(air_zoo), "obs,",
    round(diff(range(index(air_zoo))) / 86400), "days\n")
cat("  Well 20 (diel):    ", length(well20_zoo), "obs,",
    round(diff(range(index(well20_zoo))) / 86400), "days\n")
cat("  Well 21 (175m):    ", length(well21_zoo), "obs,",
    round(diff(range(index(well21_zoo))) / 86400), "days\n")
cat("  Well 10 (955m):    ", length(well10_zoo), "obs,",
    round(diff(range(index(well10_zoo))) / 86400), "days\n")
cat("  Spring 6:          ", length(spring6_zoo), "obs,",
    round(diff(range(index(spring6_zoo))) / 86400), "days\n")
cat("\n  Full individual records:\n")
cat("  Channel 3 (full):", length(stream_full_zoo), "obs,",
    round(diff(range(index(stream_full_zoo))) / 86400), "days\n")
cat("  Well 10 (full):  ", length(well10_full_zoo), "obs,",
    round(diff(range(index(well10_full_zoo))) / 86400), "days\n")
cat("  Spring 6 (full): ", length(spring6_full_zoo), "obs,",
    round(diff(range(index(spring6_full_zoo))) / 86400), "days\n")


# === 2. Overview: raw time series =============================================
#
# Visual survey of all six signals. The gradient from stream → well 20 →
# well 21 → well 10 → spring shows progressive attenuation of the diel
# cycle and damping of the seasonal amplitude.

try({
par(mfrow = c(3, 2), mar = c(3, 3, 2, 0.5), mgp = c(2, 0.7, 0))

plot_zoo_days <- function(z, main, col = "steelblue") {
  t_days <- (index(z) - index(z)[1]) / 86400
  plot(t_days, coredata(z), type = "l", col = col,
       xlab = "Days from start", ylab = "Temp (\u00B0C)", main = main)
}

plot_zoo_days(stream_zoo,  "Channel 3 (stream, 0 m)")
plot_zoo_days(air_zoo,     "UP211 air temperature", col = "firebrick")
plot_zoo_days(well20_zoo,  "Well 20 (shallow, diel detectable)")
plot_zoo_days(well21_zoo,  "Well 21 (175 m flowpath)")
plot_zoo_days(well10_zoo,  "Well 10 (955 m flowpath)")
plot_zoo_days(spring6_zoo, "Offchannel spring 6", col = "darkgreen")


}, silent = TRUE)
# === 3. Power spectra =========================================================
#
# What frequencies are present in each signal? The stream and air should
# show clear diel (24h) and seasonal peaks. Well 20 should show an
# attenuated diel peak. Wells 21 and 10 should show seasonal only.
# Spring should show a weak seasonal peak at best.

compute_power <- function(z, label) {
  vals <- coredata(z)
  vals[is.na(vals)] <- mean(vals, na.rm = TRUE)  # fill NAs for FFT
  n <- length(vals)
  dt <- median(diff(index(z)))
  fft_res <- fft(vals - mean(vals))
  freqs <- (0:(n - 1)) / (n * dt)
  power <- Mod(fft_res)^2 / n
  periods_hr <- 1 / freqs / 3600
  nyq <- floor(n / 2)
  # Subset to positive frequencies, drop DC and Nyquist edge
  p_hr <- periods_hr[2:nyq]
  pwr  <- power[2:nyq]
  # Guard against Inf/NaN from zero-frequency leakage
  keep <- is.finite(p_hr) & is.finite(pwr) & p_hr > 0 & pwr > 0
  list(periods_hr = p_hr[keep], power = pwr[keep], label = label)
}

spectra <- list(
  compute_power(stream_full_zoo, "Stream (907d)"),
  compute_power(air_zoo, "Air"),
  compute_power(well20_zoo, "Well 20"),
  compute_power(well21_zoo, "Well 21"),
  compute_power(well10_full_zoo, "Well 10 (6yr)"),
  compute_power(spring6_full_zoo, "Spring 6 (449d)")
)

try({
par(mfrow = c(3, 2), mar = c(3, 3, 2, 0.5), mgp = c(2, 0.7, 0))
for (sp in spectra) {
  plot(sp$periods_hr, sp$power, type = "l", log = "xy",
       xlab = "Period (hours)", ylab = "Power",
       main = paste("Power spectrum \u2014", sp$label), col = "steelblue")
  abline(v = 24, col = "red", lty = 2)
  abline(v = 12, col = "orange", lty = 2)
  text(24, max(sp$power) * 0.1, "24h", col = "red", pos = 4, cex = 0.7)
  text(12, max(sp$power) * 0.03, "12h", col = "orange", pos = 4, cex = 0.7)
}


}, silent = TRUE)
# === 4. Channel 3: diel fits ==================================================
#
# The stream has the strongest diel signal. Detrend the seasonal cycle,
# window into 7-day chunks, and fit each with OLS and FFT.

# Use continuous concurrent record for windowed diel analysis (full record has 394-day gap)
stream_dt <- detrend_temperature(stream_zoo, method = "loess", span = 0.05,
                                  return_trend = TRUE)

# Wrap in single-column zoo with name (fit_ols/fit_fft need column names)
name_zoo <- function(z, nm) {
  zoo(matrix(coredata(z), ncol = 1, dimnames = list(NULL, nm)),
      order.by = index(z))
}

stream_wins <- stream_dt$detrended %>%
  name_zoo("stream") %>%
  window_temperature(width = 7 * 86400, step = 86400)

stream_win_fits <- stream_wins %>%
  mutate(
    ols = map(data, ~ fit_ols(.x, 0, 86400, c(-1/8, 7/8), 10, 1)),
    fft = map(data, ~ fit_fft(.x, 0, 86400, c(-1/8, 7/8), 10, 1)),
    amp_ols   = map_dbl(ols, ~ attr(.x, "amplitudes")[1]),
    amp_fft   = map_dbl(fft, ~ attr(.x, "amplitudes")[1]),
    phase_ols = map_dbl(ols, ~ attr(.x, "phases")[1]),
    phase_fft = map_dbl(fft, ~ attr(.x, "phases")[1])
  )

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
center_days <- stream_win_fits$window_center / 86400

plot(center_days, stream_win_fits$amp_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "Days from start", ylab = "Diel amplitude (\u00B0C)",
     main = "Channel 3 \u2014 time-varying diel amplitude",
     ylim = range(c(stream_win_fits$amp_ols, stream_win_fits$amp_fft), na.rm = TRUE))
lines(center_days, stream_win_fits$amp_fft, col = "firebrick", lwd = 2, lty = 2)
legend("topright", c("OLS", "FFT"), col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = 2, cex = 0.7, bty = "n")

plot(center_days, stream_win_fits$phase_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "Days from start", ylab = "Diel phase (radians)",
     main = "Channel 3 \u2014 time-varying diel phase")
lines(center_days, stream_win_fits$phase_fft, col = "firebrick", lwd = 2, lty = 2)
legend("topright", c("OLS", "FFT"), col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = 2, cex = 0.7, bty = "n")

cat("\n=== Channel 3 diel amplitude ===\n")
cat("  Summer peak:", round(max(stream_win_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")
cat("  Winter trough:", round(min(stream_win_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")


# === 5. Channel 3: annual fit =================================================
#
# Fit the annual cycle on the full record (daily averages). OLS vs FFT
# on a ~15-month record — the annual bin won't land exactly at 365.25 days.

stream_daily <- stream_full_zoo %>%
  {zoo(coredata(.), order.by = index(.))} %>%
  # Aggregate to daily means
  aggregate(by = as.numeric(index(.)) %/% 86400 * 86400, FUN = mean, na.rm = TRUE)

stream_daily_named <- name_zoo(stream_daily, "stream")

ols_annual <- fit_ols(stream_daily_named, mean(coredata(stream_daily)),
                       365.25 * 86400, c(-1/8, 7/8), 10, 1)
fft_annual <- fit_fft(stream_daily_named, mean(coredata(stream_daily)),
                       365.25 * 86400, c(-1/8, 7/8), 10, 1)

cat("\n=== Channel 3 annual fit ===\n")
cat("  OLS amplitude:", round(attr(ols_annual, "amplitudes"), 2), "\u00B0C\n")
cat("  FFT amplitude:", round(attr(fft_annual, "amplitudes"), 2), "\u00B0C\n")
cat("  OLS phase:", round(attr(ols_annual, "phases"), 3), "rad\n")
cat("  FFT phase:", round(attr(fft_annual, "phases"), 3), "rad\n")


# === 6. Air temperature: diel + seasonal ======================================
#
# Air has bigger diel swings than the stream but also synoptic noise
# (weather fronts). How does the time-varying diel amplitude compare
# to the stream? The stream integrates (filters) the air signal.

air_dt <- detrend_temperature(air_zoo, method = "loess", span = 0.05)

air_wins <- air_dt %>%
  name_zoo("air") %>%
  window_temperature(width = 7 * 86400, step = 86400)

air_win_fits <- air_wins %>%
  mutate(
    amp_ols = map_dbl(data, function(d) {
      fit <- fit_ols(d, 0, 86400, c(-1/8, 7/8), 10, 1)
      attr(fit, "amplitudes")[1]
    })
  )

# Compare stream and air diel amplitude on the same plot
# Align by date (air starts ~3 weeks later)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

stream_dates <- as.POSIXct(stream_win_fits$window_center, origin = "1970-01-01", tz = "UTC")
air_dates    <- as.POSIXct(air_win_fits$window_center, origin = "1970-01-01", tz = "UTC")

plot(stream_dates, stream_win_fits$amp_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "", ylab = "Diel amplitude (\u00B0C)",
     main = "Diel amplitude: stream vs air",
     ylim = c(0, max(air_win_fits$amp_ols, na.rm = TRUE) * 1.1))
lines(air_dates, air_win_fits$amp_ols, col = "firebrick", lwd = 2)
legend("topright", c("Stream (channel 3)", "Air (UP211)"),
       col = c("steelblue", "firebrick"), lwd = 2, cex = 0.7, bty = "n")

cat("\n=== Air vs Stream diel amplitude ===\n")
cat("  Air summer peak:", round(max(air_win_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")
cat("  Stream summer peak:", round(max(stream_win_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")
cat("  The stream damps the air signal substantially.\n")


# === 7. Well 20: diel detection at depth ======================================
#
# Well 20 is a shallow floodplain well (~63m from channel 3) with a
# detectable diel signal in summer. This section explores the limits
# of diel detection: when does the fitted amplitude become noise?

well20_dt <- detrend_temperature(well20_zoo, method = "loess", span = 0.1)

well20_wins <- well20_dt %>%
  name_zoo("well20") %>%
  window_temperature(width = 7 * 86400, step = 86400)

well20_win_fits <- well20_wins %>%
  mutate(
    amp_ols = map_dbl(data, function(d) {
      fit <- fit_ols(d, 0, 86400, c(-1/8, 7/8), 10, 1)
      attr(fit, "amplitudes")[1]
    }),
    phase_ols = map_dbl(data, function(d) {
      fit <- fit_ols(d, 0, 86400, c(-1/8, 7/8), 10, 1)
      attr(fit, "phases")[1]
    })
  )

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
w20_dates <- as.POSIXct(well20_win_fits$window_center, origin = "1970-01-01", tz = "UTC")

plot(w20_dates, well20_win_fits$amp_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "", ylab = "Diel amplitude (\u00B0C)",
     main = "Well 20 \u2014 diel amplitude (\u2248 63 m from channel)")
abline(h = 0, lty = 2, col = "gray50")

plot(w20_dates, well20_win_fits$phase_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "", ylab = "Diel phase (radians)",
     main = "Well 20 \u2014 diel phase")

cat("\n=== Well 20 diel amplitude ===\n")
cat("  Summer peak:", round(max(well20_win_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")
cat("  Winter trough:", round(min(well20_win_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")
cat("  The diel signal emerges in summer and vanishes by late fall.\n")
cat("  When amplitude is near zero, the phase becomes meaningless (noise).\n")


# === 8. Wells 21 and 10: bashing against the diel limit =======================
#
# These wells are at 175m and 955m flowpath distance. The diel signal
# should be fully attenuated. We fit diel anyway — the point is to show
# what "no signal" looks like through the fitting functions.
# "Zeroes are data too."

diel_limit_test <- function(z, label) {
  dt <- detrend_temperature(z, method = "loess", span = 0.1)
  wins <- dt %>%
    name_zoo(label) %>%
    window_temperature(width = 7 * 86400, step = 86400)
  wins %>%
    mutate(
      amp_ols = map_dbl(data, function(d) {
        fit <- fit_ols(d, 0, 86400, c(-1/8, 7/8), 10, 1)
        attr(fit, "amplitudes")[1]
      })
    )
}

well21_diel <- diel_limit_test(well21_zoo, "well21")
well10_diel <- diel_limit_test(well10_zoo, "well10")

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

w21_dates <- as.POSIXct(well21_diel$window_center, origin = "1970-01-01", tz = "UTC")
w10_dates <- as.POSIXct(well10_diel$window_center, origin = "1970-01-01", tz = "UTC")

# Same y-axis as well 20 for direct comparison
y_max <- max(well20_win_fits$amp_ols, na.rm = TRUE) * 1.1

plot(w20_dates, well20_win_fits$amp_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "", ylab = "Diel amplitude (\u00B0C)",
     main = "Diel amplitude: well 20 vs well 21 vs well 10",
     ylim = c(0, y_max),
     xlim = range(c(w20_dates, w21_dates, w10_dates)))
lines(w21_dates, well21_diel$amp_ols, col = "orange", lwd = 2)
lines(w10_dates, well10_diel$amp_ols, col = "red", lwd = 2)
legend("topright", c("Well 20 (shallow)", "Well 21 (175m)", "Well 10 (955m)"),
       col = c("steelblue", "orange", "red"), lwd = 2, cex = 0.7, bty = "n")

# Zoom in on the noise floor
plot(w21_dates, well21_diel$amp_ols, type = "l", col = "orange", lwd = 2,
     xlab = "", ylab = "Diel amplitude (\u00B0C)",
     main = "Diel amplitude: noise floor (wells 21 + 10)",
     ylim = c(0, max(c(well21_diel$amp_ols, well10_diel$amp_ols), na.rm = TRUE) * 1.1))
lines(w10_dates, well10_diel$amp_ols, col = "red", lwd = 2)
legend("topright", c("Well 21 (175m)", "Well 10 (955m)"),
       col = c("orange", "red"), lwd = 2, cex = 0.7, bty = "n")

cat("\n=== Diel detection limits ===\n")
cat("  Well 20 max diel amp:", round(max(well20_win_fits$amp_ols, na.rm = TRUE), 3), "\u00B0C\n")
cat("  Well 21 max diel amp:", round(max(well21_diel$amp_ols, na.rm = TRUE), 3), "\u00B0C\n")
cat("  Well 10 max diel amp:", round(max(well10_diel$amp_ols, na.rm = TRUE), 3), "\u00B0C\n")
cat("  Well 21 and 10 amplitudes are at sensor resolution / noise floor.\n")


# === 9. Spring 6: seasonal fit ================================================
#
# The spring represents the deep groundwater end-member. No diel signal
# (sub-resolution). Using the full 449-day record (> 1 year) for a more
# reliable annual fit than the original 248-day excerpt.

spring6_full_named <- name_zoo(spring6_full_zoo, "spring6")

ols_spring_annual <- fit_ols(spring6_full_named, mean(coredata(spring6_full_zoo), na.rm = TRUE),
                              365.25 * 86400, c(-1/8, 7/8), 10, 1)
fft_spring_annual <- fit_fft(spring6_full_named, mean(coredata(spring6_full_zoo), na.rm = TRUE),
                              365.25 * 86400, c(-1/8, 7/8), 10, 1)

# Also try diel — should be near zero
ols_spring_diel <- fit_ols(spring6_full_named, mean(coredata(spring6_full_zoo), na.rm = TRUE),
                            86400, c(-1/8, 7/8), 10, 1)

cat("\n=== Spring 6 fits (full 449-day record) ===\n")
cat("  Annual amplitude (OLS):", round(attr(ols_spring_annual, "amplitudes"), 3), "\u00B0C\n")
cat("  Annual amplitude (FFT):", round(attr(fft_spring_annual, "amplitudes"), 3), "\u00B0C\n")
cat("  Annual phase (OLS):", round(attr(ols_spring_annual, "phases"), 3), "rad\n")
cat("  Diel amplitude (OLS):",   round(attr(ols_spring_diel, "amplitudes"), 4), "\u00B0C\n")
cat("  This is the thermal signature of deep groundwater discharge.\n")

# Visualize: spring raw data with fit_ols annual overlay
# The fitted values from fit_ols are the OLS cosine reconstruction
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

# Reconstruct the OLS cosine fit from amplitude + phase
t_sec <- index(spring6_full_zoo) - index(spring6_full_zoo)[1]
omega_annual <- 2 * pi / (365.25 * 86400)
amp_s  <- attr(ols_spring_annual, "amplitudes")[1]
pha_s  <- attr(ols_spring_annual, "phases")[1]
mean_s <- mean(coredata(spring6_full_zoo), na.rm = TRUE)
fitted_annual <- mean_s + amp_s * cos(omega_annual * t_sec - pha_s)

spring_dates <- as.POSIXct(index(spring6_full_zoo), origin = "1970-01-01", tz = "UTC")
plot(spring_dates, coredata(spring6_full_zoo), type = "l", col = "darkgreen",
     xlab = "", ylab = "Temp (\u00B0C)",
     main = "Offchannel spring 6 \u2014 fit_ols annual overlay (449 days)")
lines(spring_dates, fitted_annual, col = "red", lwd = 2)
legend("topleft", c("Observed", "OLS annual fit"),
       col = c("darkgreen", "red"), lty = 1, lwd = c(1, 2), cex = 0.7, bty = "n")

# === 10. Annual inverse solution: 3-site analysis =============================
#
# This is the full thObservedSeries() call on daily-averaged data from
# channel 3 (0m), well 21 (175m), and well 10 (955m). Annual period.
# Recovers eta, Darcy flux, thermal diffusivity, dispersivity.
#
# Flowpath lengths from field mapping:
#   channel 3 → well 21: 175 m
#   channel 3 → well 10: 955 m

xVals <- c("channel_3_0m" = 0, "well_21_175m" = 175, "well_10_955m" = 955)

# Build multi-column daily zoo
three_site_zoo <- suppressWarnings(zoo(
  as.matrix(three_site[, c("channel_3_0m", "well_21_175m", "well_10_955m")]),
  order.by = as.numeric(three_site$date_time)
))

# Aggregate to daily means
daily_idx <- as.numeric(index(three_site_zoo)) %/% 86400 * 86400
three_site_daily <- aggregate(three_site_zoo, by = daily_idx, FUN = mean, na.rm = TRUE)

# Interpolate remaining NAs
three_site_daily <- na.approx(three_site_daily, na.rm = FALSE)

cat("\n=== 3-site daily data for inverse model ===\n")
cat("  Rows:", nrow(three_site_daily), "\n")
cat("  Date range:", as.character(as.POSIXct(range(index(three_site_daily)),
    origin = "1970-01-01", tz = "UTC")), "\n")

# Aquifer properties (basalt defaults from temperheic)
aquifer <- generate_example_aquifer()

# Run inverse model: OLS
inv_ols <- thObservedSeries(
  empiricalData = three_site_daily,
  xVals = xVals,
  aquifer = aquifer,
  boundaryMean = mean(coredata(three_site_daily[, "channel_3_0m"]), na.rm = TRUE),
  period = 365.25 * 86400,
  headGrad = 0.008,
  nmin = 180,
  fit_method = "ols"
)

# Run inverse model: FFT
inv_fft <- thObservedSeries(
  empiricalData = three_site_daily,
  xVals = xVals,
  aquifer = aquifer,
  boundaryMean = mean(coredata(three_site_daily[, "channel_3_0m"]), na.rm = TRUE),
  period = 365.25 * 86400,
  headGrad = 0.008,
  nmin = 180,
  fit_method = "fft"
)

cat("\n=== Annual inverse results (OLS) ===\n")
cat("\nAmplitudes (\u00B0C):\n")
print(round(inv_ols$amplitude, 3))
cat("\nAmplitude ratios:\n")
print(round(inv_ols$ampRatio, 4))
cat("\nPhase differences (radians):\n")
print(round(inv_ols$deltaPhaseRadians, 4))
cat("\nEta:\n")
print(round(inv_ols$eta, 4))
cat("\nDarcy flux (m/s):\n")
print(inv_ols$darcyFlux)
cat("\nEffective diffusivity (m2/s):\n")
print(inv_ols$diffusivity_effective_empirical)
cat("\nDispersivity (m):\n")
print(inv_ols$dispersivity)
cat("\nHydraulic conductivity (m/s):\n")
print(inv_ols$hydraulicCond)

cat("\n=== OLS vs FFT comparison ===\n")
cat("Amplitude ratio [1,2] — OLS:", round(inv_ols$ampRatio[1, 2], 4),
    "FFT:", round(inv_fft$ampRatio[1, 2], 4), "\n")
cat("Delta phase [1,2] — OLS:", round(inv_ols$deltaPhaseRadians[1, 2], 4),
    "FFT:", round(inv_fft$deltaPhaseRadians[1, 2], 4), "\n")
cat("Darcy flux [1,2] — OLS:", format(inv_ols$darcyFlux[1, 2], digits = 4),
    "FFT:", format(inv_fft$darcyFlux[1, 2], digits = 4), "\n")

# Plot: observed daily temperatures with annual cosine fits
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
ts_dates <- as.POSIXct(index(three_site_daily), origin = "1970-01-01", tz = "UTC")

plot(ts_dates, coredata(three_site_daily[, 1]), type = "l", col = "steelblue",
     ylim = range(three_site_daily, na.rm = TRUE),
     xlab = "", ylab = "Daily mean temp (\u00B0C)",
     main = "3-site annual temperature \u2014 observed")
lines(ts_dates, coredata(three_site_daily[, 2]), col = "orange")
lines(ts_dates, coredata(three_site_daily[, 3]), col = "red")
legend("topright", c("Channel 3 (0m)", "Well 21 (175m)", "Well 10 (955m)"),
       col = c("steelblue", "orange", "red"), lty = 1, cex = 0.7, bty = "n")


cat("\n=== Script complete. Explore: stream_zoo, air_zoo, well20_zoo,\n")
cat("    well21_zoo, well10_zoo, spring6_zoo, stream_win_fits,\n")
cat("    air_win_fits, well20_win_fits, well21_diel, well10_diel,\n")
cat("    inv_ols, inv_fft ===\n")

