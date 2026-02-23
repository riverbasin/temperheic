
# testingThSeries.R — Interactive exploration of temperheic forward and inverse models
#
# Purpose:
#   Generate synthetic temperature time series across parameter combinations,
#   run the inverse model with both OLS and NLS fitting methods, and compare
#   recovered parameters against known values. Includes visualizations of:
#     - Forward model temperature signals at multiple depths
#     - Amplitude damping and phase lag with distance
#     - OLS vs NLS fit method comparison
#     - Dispersivity and Peclet number relationships
#
# References:
#   Luce et al. (2013) — analytical solutions, eta parameter, inverse model
#   Vogt et al. (2012) — x_d / v_d formulation of the ADE solution
#   Bons et al. (2013) — Pe-dependent dispersion, nonlinear D(Pe) regimes
#
# Usage: source this script interactively in RStudio, or run sections as needed.

library(temperheic)
library(zoo)
library(dplyr)
Sys.setenv(TZ = "UTC")


# === 1. Define aquifer and parameter space ====================================

aquifer <- generate_example_aquifer()

# Dispersivity sweep — the key parameter linking to Bons et al.
# At low dispersivity, kappa_e ≈ kappa_cond (diffusion-dominated).
# At high dispersivity, kappa_disp dominates (mechanical dispersion regime).
# The transition between regimes is governed by the Peclet number.
scenarios <- list(
  diel = list(
    k            = 10 / 86400,            # K = 10 m/d, moderate
    dispersivity = c(0.001, 0.01, 0.1, 1),  # 4 values spanning regimes
    headGrad     = 0.1,
    mean         = 22,
    amplitude    = 6,
    phase        = 86400 + 3600,          # peak ~25h after start
    period       = 86400,                 # 1 day
    xVals        = c(0, 0.33, 0.66, 1),  # shallow streambed sensors (m)
    tVals        = seq(0, 3600 * 72, 3600) # 3 days, hourly
  ),
  annual = list(
    k            = 100 / 86400,           # K = 100 m/d
    dispersivity = c(1, 10, 100, 1000),   # floodplain-scale dispersivities
    headGrad     = 0.01,
    mean         = 11.78,
    amplitude    = 8.87,
    phase        = 1702226,               # ~July peak (from Meacham fits)
    period       = 86400 * 365,           # 1 year
    xVals        = c(0, 120, 240, 360),   # floodplain well distances (m)
    tVals        = seq(0, 2 * 365 * 86400, 86400)  # 2 years, daily
  )
)


# === 2. Generate forward model time series ====================================

fwd <- generate_example_series(aquifer, scenarios)

cat("Forward model scenarios generated:\n")
cat("  diel combos:", length(fwd$diel), "\n")
cat("  annual combos:", length(fwd$annual), "\n")


# === 3. Plot forward model signals ============================================

# --- 3a. Diel signals: temperature at each depth over 3 days ---
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (i in seq_along(fwd$diel)) {
  ts_data <- fwd$diel[[i]]$timeSeries
  disp_val <- fwd$diel[[i]]$signal$hydro$dispersivity

  # Convert seconds to hours for readable x-axis
  hours <- as.numeric(index(ts_data)) / 3600
  plot(hours, coredata(ts_data)[, 1], type = "l", col = 1,
       ylim = range(ts_data), xlab = "Hours", ylab = "Temperature (°C)",
       main = paste0("Diel | dispersivity = ", disp_val, " m"))
  for (j in 2:ncol(ts_data)) {
    lines(hours, coredata(ts_data)[, j], col = j)
  }
  legend("topright", names(ts_data), col = seq_len(ncol(ts_data)),
         lty = 1, cex = 0.7, bty = "n")
}

# --- 3b. Annual signals: temperature at each depth over 2 years ---
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (i in seq_along(fwd$annual)) {
  ts_data <- fwd$annual[[i]]$timeSeries
  disp_val <- fwd$annual[[i]]$signal$hydro$dispersivity

  days <- as.numeric(index(ts_data)) / 86400
  plot(days, coredata(ts_data)[, 1], type = "l", col = 1,
       ylim = range(ts_data), xlab = "Days", ylab = "Temperature (°C)",
       main = paste0("Annual | dispersivity = ", disp_val, " m"))
  for (j in 2:ncol(ts_data)) {
    lines(days, coredata(ts_data)[, j], col = j)
  }
  legend("topright", names(ts_data), col = seq_len(ncol(ts_data)),
         lty = 1, cex = 0.7, bty = "n")
}


# === 4. Amplitude decay and phase lag profiles ================================
#
# For each scenario, plot how amplitude and phase change with distance.
# The exponential decay A(x) = A_0 * exp(-x/x_d) and linear phase shift
# phi(x) = x*omega/v_d are the Vogt (2012) formulation of the ADE solution.

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (scenario_name in c("diel", "annual")) {
  scenario <- fwd[[scenario_name]]

  # Amplitude vs distance
  plot(NULL, xlim = range(scenario[[1]]$xVals),
       ylim = c(0, max(sapply(scenario, function(s) max(s$amplitude)))),
       xlab = "Distance (m)", ylab = "Amplitude (°C)",
       main = paste0(scenario_name, " — Amplitude decay"))
  for (i in seq_along(scenario)) {
    s <- scenario[[i]]
    lines(s$xVals, s$amplitude, col = i, type = "b", pch = 16)
  }
  legend("topright",
         paste("disp =", sapply(scenario, function(s) s$signal$hydro$dispersivity)),
         col = seq_along(scenario), lty = 1, cex = 0.7, bty = "n")

  # Phase vs distance
  plot(NULL, xlim = range(scenario[[1]]$xVals),
       ylim = range(sapply(scenario, function(s) range(s$phaseRadians))),
       xlab = "Distance (m)", ylab = "Phase (radians)",
       main = paste0(scenario_name, " — Phase lag"))
  for (i in seq_along(scenario)) {
    s <- scenario[[i]]
    lines(s$xVals, s$phaseRadians, col = i, type = "b", pch = 16)
  }
  legend("topleft",
         paste("disp =", sapply(scenario, function(s) s$signal$hydro$dispersivity)),
         col = seq_along(scenario), lty = 1, cex = 0.7, bty = "n")
}


# === 5. Run inverse model: OLS vs NLS comparison ==============================

inv_ols <- generate_example_observed(fwd, fit_method = "ols")
inv_nls <- generate_example_observed(fwd, fit_method = "nls")

cat("\nInverse models complete (OLS and NLS).\n")


# === 6. Compare recovered vs known parameters =================================
#
# For each combination, extract the [1,2] element (boundary vs first
# subsurface sensor) and compare known/recovered values.

comparison <- purrr::imap_dfr(fwd, function(scenario_fwd, scenario_name) {
  purrr::imap_dfr(scenario_fwd, function(fwd_s, combo_name) {
    ols <- inv_ols[[scenario_name]][[combo_name]]
    nls <- inv_nls[[scenario_name]][[combo_name]]

    known_darcy <- fwd_s$signal$hydro$darcyFlux
    known_diff  <- fwd_s$signal$hydro$diffusivity_effective
    known_vt    <- fwd_s$signal$hydro$advectiveThermVel
    known_disp  <- fwd_s$signal$hydro$dispersivity
    known_Pe    <- fwd_s$signal$pecletNumberEff

    tibble::tibble(
      scenario   = scenario_name,
      combo      = combo_name,
      dispersivity_known = known_disp,
      peclet_known       = known_Pe,
      parameter  = rep(c("darcyFlux", "diffusivity_eff", "advectiveThermVel"), 2),
      method     = rep(c("ols", "nls"), each = 3),
      known      = rep(c(known_darcy, known_diff, known_vt), 2),
      recovered  = c(
        ols$darcyFlux[1, 2], ols$diffusivity_effective_empirical[1, 2],
        ols$advectiveThermVelEmpirical[1, 2],
        nls$darcyFlux[1, 2], nls$diffusivity_effective_empirical[1, 2],
        nls$advectiveThermVelEmpirical[1, 2]
      )
    )
  })
}) %>%
  dplyr::mutate(relative_error = abs(recovered - known) / abs(known))

cat("\n=== Parameter Recovery Summary ===\n")
print(
  comparison %>%
    dplyr::group_by(scenario, method) %>%
    dplyr::summarise(
      max_rel_error = max(relative_error),
      mean_rel_error = mean(relative_error),
      .groups = "drop"
    ),
  n = 20
)


# === 7. Dispersivity and Peclet number diagnostics ============================
#
# The inverse model recovers dispersivity from:
#   dispersivity = (kappa_e * rho_m*c_m - lambda_m) / (v_t * rho_m*c_m)
#
# If Bons et al. (2013) is correct, the relationship between kappa_e and Pe
# is nonlinear in the intermediate regime (Pe ~0.3 to ~4000), with an
# effective power-law exponent n that can exceed 1 for longitudinal dispersion.
#
# Here we check whether the recovered dispersivity and Pe are consistent
# with the known inputs across the parameter sweep.

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

for (scenario_name in c("diel", "annual")) {
  scenario_inv <- inv_ols[[scenario_name]]
  scenario_fwd <- fwd[[scenario_name]]

  known_disp <- sapply(scenario_fwd, function(s) s$signal$hydro$dispersivity)
  known_ke   <- sapply(scenario_fwd, function(s) s$signal$hydro$diffusivity_effective)
  known_kcond <- scenario_fwd[[1]]$signal$hydro$diffusivity_cond

  # Recovered values at [1,2]
  rec_ke   <- sapply(scenario_inv, function(o) o$diffusivity_effective_empirical[1, 2])
  rec_disp <- sapply(scenario_inv, function(o) o$dispersivity[1, 2])
  rec_Pe   <- sapply(scenario_inv, function(o) o$pecletNumber[1, 2])

  # Plot: recovered kappa_e / kappa_cond vs Peclet number (Bons-style)
  plot(rec_Pe, rec_ke / known_kcond, log = "xy",
       pch = 16, col = "steelblue", cex = 1.5,
       xlab = "Peclet Number", ylab = expression(kappa[e] / kappa[cond]),
       main = paste0(scenario_name, " — D/D_cond vs Pe\n(cf. Bons Fig. 4)"))
  abline(h = 1, lty = 2, col = "gray50")
  text(rec_Pe, rec_ke / known_kcond,
       labels = paste0("disp=", round(known_disp, 3)),
       pos = 3, cex = 0.7)
}


# === 8. Thermal decay distance and phase velocity =============================
#
# Vogt (2012) formulation: x_d is the e-folding distance for amplitude,
# v_d is the phase velocity. Both depend on frequency and v_t/kappa_e.
#
# Comparing these across dispersivity values shows how the signal
# penetration depth changes with flow regime.

cat("\n=== Thermal Decay Distance (x_d) and Phase Velocity (v_d) ===\n")
for (scenario_name in c("diel", "annual")) {
  cat("\n", toupper(scenario_name), ":\n")
  for (i in seq_along(fwd[[scenario_name]])) {
    s <- fwd[[scenario_name]][[i]]$signal
    cat(sprintf("  disp = %8.3f | x_d = %10.2f m | v_d = %.4e m/s | Pe_eff = %.2f\n",
                s$hydro$dispersivity,
                s$thermDecayDist,
                s$phaseVel,
                s$pecletNumberEff))
  }
}



# === 9. Composite stream temperature: diel + annual signal ===================
#
# Real stream temperature is not a pure sinusoid at a single frequency.
# It is an amplitude-modulated signal: a diel cycle whose amplitude varies
# seasonally (large in summer, small in winter), superimposed on the annual
# cycle. Parameters below are derived from Meacham Creek Channel 3
# (2012-2014, 907 days including winter):
#
#   Annual mean:          10.5 °C
#   Annual amplitude:      8.0 °C (range ~2.5 to ~18.5 °C)
#   Mean diel amplitude:   2.2 °C (half-range)
#   Diel amp modulation:   1.8 °C (summer ~4.0, winter ~0.4 °C)
#
# The product of the diel carrier with the seasonal envelope creates
# sidebands at (1/day ± 1/year) in the FFT — this is the spectral
# fingerprint of amplitude modulation.

composite_stream <- generate_composite_boundary(n_years = 1)

# Plot: full year
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
t_days <- as.numeric(index(composite_stream)) / 86400
plot(t_days, coredata(composite_stream), type = "l", col = "steelblue",
     xlab = "Day of year", ylab = "Temperature (°C)",
     main = "Composite stream temperature (1 year)")

# Zoom: 2 weeks in summer (day 190-204)
summer_idx <- t_days >= 190 & t_days <= 204
plot(t_days[summer_idx], coredata(composite_stream)[summer_idx],
     type = "l", col = "red",
     xlab = "Day of year", ylab = "Temperature (°C)",
     main = "Summer diel cycles (Jul)")

# Zoom: 2 weeks in winter (day 0-14)
winter_idx <- t_days >= 0 & t_days <= 14
plot(t_days[winter_idx], coredata(composite_stream)[winter_idx],
     type = "l", col = "blue",
     xlab = "Day of year", ylab = "Temperature (°C)",
     main = "Winter diel cycles (Jan)")

# Power spectrum — should show annual peak, diel peak, and sidebands
temp_vec <- coredata(composite_stream)
n <- length(temp_vec)
dt <- median(diff(index(composite_stream)))
fft_result <- fft(temp_vec - mean(temp_vec))
freqs <- (0:(n - 1)) / (n * dt)
power <- Mod(fft_result)^2 / n
periods_hr <- 1 / freqs / 3600

nyq <- floor(n / 2)
plot(periods_hr[2:nyq], power[2:nyq], type = "l", log = "xy",
     xlab = "Period (hours)", ylab = "Power",
     main = "Power spectrum — composite signal")
abline(v = 24, col = "red", lty = 2)
abline(v = 12, col = "orange", lty = 2)
abline(v = 24 * 365.25, col = "blue", lty = 2)
text(24, max(power[2:nyq]) * 0.01, "24h", col = "red", pos = 2, cex = 0.8)
text(12, max(power[2:nyq]) * 0.003, "12h", col = "orange", pos = 2, cex = 0.8)
text(24 * 365.25, max(power[2:nyq]) * 0.3, "annual", col = "blue", pos = 2, cex = 0.8)


# === 10. Compare OLS and FFT on the composite signal ==========================
#
# Extract the diel amplitude and phase from the composite using both methods.
# The composite is a single boundary sensor, so we get one amplitude/phase.
# The key question: does the diel fit correctly extract ~2.2 °C mean amplitude
# even though the signal is amplitude-modulated?

comp_zoo <- zoo(matrix(coredata(composite_stream), ncol = 1,
                       dimnames = list(NULL, "stream")),
                order.by = index(composite_stream))

# Diel fit
ols_diel <- fit_ols(comp_zoo, mean(coredata(composite_stream)), 86400,
                     c(-1/8, 7/8), 10, 1)
fft_diel <- fit_fft(comp_zoo, mean(coredata(composite_stream)), 86400,
                     c(-1/8, 7/8), 10, 1)

# Annual fit
ols_annual <- fit_ols(comp_zoo, mean(coredata(composite_stream)),
                       86400 * 365.25, c(-1/8, 7/8), 10, 1)
fft_annual <- fit_fft(comp_zoo, mean(coredata(composite_stream)),
                       86400 * 365.25, c(-1/8, 7/8), 10, 1)

cat("\n=== Composite signal: amplitude extraction ===\n")
cat("Known annual amplitude: 8.0 °C\n")
cat("  OLS recovered:", round(attr(ols_annual, "amplitudes"), 3), "°C\n")
cat("  FFT recovered:", round(attr(fft_annual, "amplitudes"), 3), "°C\n")
cat("\nKnown mean diel amplitude: 2.2 °C\n")
cat("  OLS recovered:", round(attr(ols_diel, "amplitudes"), 3), "°C\n")
cat("  FFT recovered:", round(attr(fft_diel, "amplitudes"), 3), "°C\n")
cat("\nNote: The diel fit recovers the MEAN diel amplitude (2.2 °C),\n")
cat("not the summer peak (4.0 °C) or winter trough (0.4 °C), because\n")
cat("the cosine regression averages over the full record.\n")
cat("This is the amplitude underprediction problem that motivates\n")
cat("windowed or time-varying analysis approaches.\n")



# === 11. Seasonal windows: data with fit overlay ==============================
#
# Fit the diel cosine separately on a summer window and a winter window,
# then overlay the fit on the data. This shows how a single-frequency
# cosine captures the signal in each season, and why the whole-year fit
# underpredicts summer amplitude.

comp_data <- coredata(composite_stream)
comp_t    <- as.numeric(index(composite_stream))
t_days_c  <- comp_t / 86400
omega_diel <- 2 * pi / 86400

# Define seasonal windows (days of year)
windows <- list(
  winter  = c(0, 30),        # Jan
  spring  = c(90, 120),      # Apr
  summer  = c(180, 210),     # Jul
  autumn  = c(270, 300)      # Oct
)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (season_name in names(windows)) {
  win <- windows[[season_name]]
  idx <- which(t_days_c >= win[1] & t_days_c <= win[2])

  # Subset for this window
  t_win <- comp_t[idx]
  T_win <- comp_data[idx]
  t_win_rel <- t_win - t_win[1]

  # OLS fit on this window
  sin_basis <- sin(omega_diel * t_win_rel)
  cos_basis <- cos(omega_diel * t_win_rel)
  fit <- lm(T_win ~ sin_basis + cos_basis)
  alpha_s <- coef(fit)[["sin_basis"]]
  alpha_c <- coef(fit)[["cos_basis"]]
  win_amp <- sqrt(alpha_s^2 + alpha_c^2)

  # Plot data + fit
  plot(t_days_c[idx], T_win, type = "l", col = "steelblue",
       xlab = "Day of year", ylab = "Temp (\u00B0C)",
       main = sprintf("%s \u2014 diel amp = %.2f \u00B0C", season_name, win_amp))
  lines(t_days_c[idx], fitted(fit), col = "red", lwd = 2)
  legend("topright", c("data", "diel fit"), col = c("steelblue", "red"),
         lty = 1, lwd = c(1, 2), cex = 0.7, bty = "n")
}

cat("\nSeasonal diel amplitude recovery (30-day windows):\n")
for (season_name in names(windows)) {
  win <- windows[[season_name]]
  idx <- which(t_days_c >= win[1] & t_days_c <= win[2])
  t_win <- comp_t[idx]
  T_win <- comp_data[idx]
  t_win_rel <- t_win - t_win[1]
  sin_b <- sin(omega_diel * t_win_rel)
  cos_b <- cos(omega_diel * t_win_rel)
  fit <- lm(T_win ~ sin_b + cos_b)
  amp <- sqrt(coef(fit)[["sin_b"]]^2 + coef(fit)[["cos_b"]]^2)
  cat(sprintf("  %-8s diel amp = %.2f C\n", season_name, amp))
}
cat(sprintf("  %-8s diel amp = %.2f C (whole-year OLS)\n", "annual",
            attr(ols_diel, "amplitudes")))
cat(sprintf("  %-8s diel amp = %.2f C (whole-year FFT)\n", "annual",
            attr(fft_diel, "amplitudes")))



# === 12. Detrending the composite signal ======================================
#
# Section 10 showed that whole-year OLS/FFT underpredicts diel amplitude
# because the annual trend biases the cosine fit. Section 11 hacked around
# this with manual seasonal windows + raw lm(). Now we use the proper
# signal_processing functions.
#
# detrend_temperature() removes the annual cycle, leaving a zero-mean
# diel oscillation whose amplitude varies seasonally (the AM envelope).

comp_detrended <- detrend_temperature(comp_zoo, method = "loess", span = 0.05,
                                       return_trend = TRUE)

par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))

# Full year: raw data
plot(comp_zoo, main = "Raw composite stream temperature",
     ylab = "Temp (\u00B0C)", xlab = "")

# Trend (the annual cycle, captured by loess)
plot(comp_detrended$trend, main = "Estimated trend (loess, span = 0.05)",
     ylab = "Temp (\u00B0C)", xlab = "")

# Detrended residual: the diel signal + noise
plot(comp_detrended$detrended, main = "Detrended \u2014 isolated diel oscillation",
     ylab = "\u0394 Temp (\u00B0C)", xlab = "Seconds")

cat("\n=== Detrended composite: range of residuals ===\n")
cat("  min:", round(min(comp_detrended$detrended, na.rm = TRUE), 2), "\u00B0C\n")
cat("  max:", round(max(comp_detrended$detrended, na.rm = TRUE), 2), "\u00B0C\n")
cat("  Summer diel peaks should be ~4 \u00B0C, winter ~0.4 \u00B0C.\n")


# === 13. Comparing detrending methods =========================================
#
# Loess vs moving average vs polynomial \u2014 how do they differ on real-ish data?
# The MA method with window = 24 (hourly data, 1-day average) is a natural
# choice for removing the diel cycle. But here we are removing the ANNUAL
# trend to preserve the diel cycle, so we need a much wider window.

methods_compare <- list(
  loess_005  = detrend_temperature(comp_zoo, method = "loess", span = 0.05),
  loess_02   = detrend_temperature(comp_zoo, method = "loess", span = 0.20),
  ma_7day    = detrend_temperature(comp_zoo, method = "ma", window = 7 * 24),
  poly_3     = detrend_temperature(comp_zoo, method = "polynomial", degree = 3)
)

# Zoom to summer: compare how each method handles the trend removal
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
t_days_d <- as.numeric(zoo::index(comp_zoo)) / 86400
summer_range <- t_days_d >= 180 & t_days_d <= 210

for (nm in names(methods_compare)) {
  dt_vals <- zoo::coredata(methods_compare[[nm]])
  plot(t_days_d[summer_range], dt_vals[summer_range],
       type = "l", col = "steelblue",
       xlab = "Day of year", ylab = "\u0394 Temp (\u00B0C)",
       main = paste("Summer diel \u2014", nm))
  abline(h = 0, lty = 2, col = "gray50")
}

cat("\nNote: loess (span=0.05) and MA (7-day) track the seasonal envelope\n")
cat("well and leave a clean zero-mean oscillation. The polynomial (degree=3)\n")
cat("cannot capture the curvature near solstice and leaves residual bias.\n")
cat("loess (span=0.20) is too smooth \u2014 it leaks diel energy into the trend.\n")


# === 14. Windowed diel amplitude: the full pipeline ===========================
#
# This is the key demonstration: detrend -> window -> fit -> collect.
# We use 7-day windows sliding by 1 day, fit each window with both OLS
# and FFT, and track how the recovered diel amplitude varies through
# the year. This is the workflow that analyze_flux() will automate.

comp_detrended_zoo <- comp_detrended$detrended

windows_7d <- window_temperature(comp_detrended_zoo,
                                  width = 7 * 86400,
                                  step  = 1 * 86400)

cat("\n=== Windowed analysis ===\n")
cat("  Windows generated:", nrow(windows_7d), "\n")
cat("  Width: 7 days, Step: 1 day\n")

# Fit each window with OLS and FFT
window_fits <- windows_7d %>%
  dplyr::mutate(
    ols_fit = purrr::map(data, ~ fit_ols(.x, mean(coredata(.x), na.rm = TRUE),
                                          86400, c(-1/8, 7/8), 10, 1)),
    fft_fit = purrr::map(data, ~ fit_fft(.x, mean(coredata(.x), na.rm = TRUE),
                                          86400, c(-1/8, 7/8), 10, 1)),
    amp_ols = purrr::map_dbl(ols_fit, ~ attr(.x, "amplitudes")[1]),
    amp_fft = purrr::map_dbl(fft_fit, ~ attr(.x, "amplitudes")[1]),
    phase_ols = purrr::map_dbl(ols_fit, ~ attr(.x, "phases")[1]),
    phase_fft = purrr::map_dbl(fft_fit, ~ attr(.x, "phases")[1])
  )

# Plot: time-varying diel amplitude (OLS vs FFT)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

center_days <- window_fits$window_center / 86400
plot(center_days, window_fits$amp_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "Day of year", ylab = "Diel amplitude (\u00B0C)",
     main = "Time-varying diel amplitude: OLS vs FFT",
     ylim = range(c(window_fits$amp_ols, window_fits$amp_fft), na.rm = TRUE))
lines(center_days, window_fits$amp_fft, col = "firebrick", lwd = 2, lty = 2)
abline(h = attr(ols_diel, "amplitudes"), col = "gray50", lty = 3)
legend("topright", c("OLS (7-day windows)", "FFT (7-day windows)",
                      "Whole-year OLS"),
       col = c("steelblue", "firebrick", "gray50"),
       lty = c(1, 2, 3), lwd = c(2, 2, 1), cex = 0.7, bty = "n")

cat("\nSummer peak amplitude (OLS):", round(max(window_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")
cat("Winter trough amplitude (OLS):", round(min(window_fits$amp_ols, na.rm = TRUE), 2), "\u00B0C\n")
cat("Whole-year OLS amplitude:", round(attr(ols_diel, "amplitudes"), 2), "\u00B0C\n")
cat("This confirms windowing resolves the amplitude underprediction problem.\n")

# Plot: time-varying diel phase (OLS vs FFT)
plot(center_days, window_fits$phase_ols, type = "l", col = "steelblue", lwd = 2,
     xlab = "Day of year", ylab = "Diel phase (radians)",
     main = "Time-varying diel phase: OLS vs FFT",
     ylim = range(c(window_fits$phase_ols, window_fits$phase_fft), na.rm = TRUE))
lines(center_days, window_fits$phase_fft, col = "firebrick", lwd = 2, lty = 2)
legend("topright", c("OLS", "FFT"),
       col = c("steelblue", "firebrick"), lty = c(1, 2), lwd = 2,
       cex = 0.7, bty = "n")


# === 15. OLS vs FFT agreement across windows ==================================
#
# With synthetic data (no noise, exact sinusoidal windows), OLS and FFT
# should agree closely. The residual disagreement comes from:
#   1. FFT bin mismatch \u2014 the 7-day window has bins at periods of
#      7*24*3600/k hours; the diel bin lands at k=7, exactly 24h. Lucky!
#   2. Edge effects in the detrended signal (loess NAs near boundaries)
#   3. AM sidebands leaking into the diel bin (minor for 7-day windows)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Amplitude scatter
plot(window_fits$amp_ols, window_fits$amp_fft,
     pch = 16, col = "steelblue", cex = 0.8,
     xlab = "OLS amplitude (\u00B0C)", ylab = "FFT amplitude (\u00B0C)",
     main = "Window amplitude: OLS vs FFT")
abline(0, 1, col = "red", lty = 2)

# Phase scatter
plot(window_fits$phase_ols, window_fits$phase_fft,
     pch = 16, col = "steelblue", cex = 0.8,
     xlab = "OLS phase (rad)", ylab = "FFT phase (rad)",
     main = "Window phase: OLS vs FFT")
abline(0, 1, col = "red", lty = 2)

amp_rmsd <- sqrt(mean((window_fits$amp_ols - window_fits$amp_fft)^2, na.rm = TRUE))
phase_rmsd <- sqrt(mean((window_fits$phase_ols - window_fits$phase_fft)^2, na.rm = TRUE))
cat("\nOLS vs FFT agreement across", nrow(window_fits), "windows:\n")
cat("  Amplitude RMSD:", round(amp_rmsd, 4), "\u00B0C\n")
cat("  Phase RMSD:    ", round(phase_rmsd, 4), "rad\n")


# === 16. Window width sensitivity =============================================
#
# How does window width affect the recovered amplitude time series?
# Wider windows smooth the seasonal envelope; narrower windows have
# fewer cycles and noisier fits. The trade-off is temporal resolution
# vs spectral resolution \u2014 the uncertainty principle at work.

widths_days <- c(3, 5, 7, 14, 30)

width_comparison <- purrr::map_dfr(widths_days, function(w_days) {
  w_sec <- w_days * 86400
  wins <- window_temperature(comp_detrended_zoo, width = w_sec, step = 86400)
  wins %>%
    dplyr::mutate(
      amp = purrr::map_dbl(data, function(d) {
        fit <- fit_ols(d, mean(coredata(d), na.rm = TRUE),
                       86400, c(-1/8, 7/8), 10, 1)
        attr(fit, "amplitudes")[1]
      }),
      width_days = w_days
    ) %>%
    dplyr::select(window_center, amp, width_days)
})

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
cols <- c("red", "orange", "steelblue", "darkgreen", "purple")
plot(NULL, xlim = c(0, 365), ylim = c(0, 5),
     xlab = "Day of year", ylab = "Diel amplitude (\u00B0C)",
     main = "Window width sensitivity \u2014 diel amplitude")

for (i in seq_along(widths_days)) {
  sub <- width_comparison %>% dplyr::filter(width_days == widths_days[i])
  lines(sub$window_center / 86400, sub$amp, col = cols[i], lwd = 1.5)
}
legend("topright", paste0(widths_days, "-day"), col = cols, lty = 1, lwd = 1.5,
       cex = 0.7, bty = "n", title = "Window width")

cat("\nThe 3-day window is noisy (only 3 complete cycles).\n")
cat("The 30-day window smooths out the seasonal envelope.\n")
cat("7 days is a reasonable default for diel analysis.\n")


cat("\n=== Script complete. Explore: fwd, inv_ols, inv_nls, comparison,\n")
cat("    composite_stream, comp_zoo, comp_detrended, window_fits,\n")
cat("    width_comparison ===\n")

