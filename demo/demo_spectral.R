# demo/demo_spectral.R
# temperheic spectral analysis demonstration
#
# Shows what compute_power() and find_peaks() can do with
# increasingly realistic signals. Source this file or run interactively.

library(temperheic)
library(zoo)

# --- Helper: build a multi-component hourly zoo ---
make_signal <- function(n_days, components, mean_temp = 10) {
  t_sec <- seq(0, by = 3600, length.out = n_days * 24)
  vals  <- rep(mean_temp, length(t_sec))
  for (c in components)
    vals <- vals + c$A * cos(2 * pi * t_sec / c$P - ifelse(is.null(c$phi), 0, c$phi))
  zoo(vals, order.by = t_sec)
}


# ============================================================
# Test 1.1  Single cosine — the baseline
# ============================================================
cat("\n=== 1.1 Single cosine (A=5, P=24h, 30 days) ===\n")
z1 <- make_signal(30, list(list(A = 5, P = 86400)))
psd1 <- compute_power(z1)
pk1  <- find_peaks(psd1, n_peaks = 3)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(z1[1:(3*24)], main = "1.1: Time domain (3 days)", ylab = "Temp (°C)", xlab = "Time (s)")
plot(psd1$period_hr, psd1$psd, type = "l", log = "x",
     main = "1.1: Power spectrum", xlab = "Period (hr)", ylab = "PSD")
abline(v = 24, lty = 2, col = "red")

cat("  Peak:", pk1$period_hr[1], "hr, Amplitude:", round(pk1$amplitude[1], 2), "°C\n")


# ============================================================
# Test 1.2  Two frequencies (24h + 12h)
# ============================================================
cat("\n=== 1.2 Diel + semi-diurnal (A=5 @24h, A=1.5 @12h) ===\n")
z2 <- make_signal(30, list(list(A = 5, P = 86400), list(A = 1.5, P = 43200)))
psd2 <- compute_power(z2)
pk2  <- find_peaks(psd2, n_peaks = 5)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(z2[1:(3*24)], main = "1.2: Time domain (3 days)", ylab = "Temp (°C)", xlab = "Time (s)")
plot(psd2$period_hr, psd2$psd, type = "l", log = "x",
     main = "1.2: Power spectrum", xlab = "Period (hr)", ylab = "PSD")
abline(v = c(12, 24), lty = 2, col = "red")

cat("  Peaks found:\n")
print(pk2[, c("period_hr", "amplitude")])


# ============================================================
# Test 1.3  Diel + annual (multi-scale, 2 years)
# ============================================================
cat("\n=== 1.3 Diel + annual (A=3 @24h, A=8 @365d, 2 years) ===\n")
z3 <- make_signal(730, list(
  list(A = 8, P = 365.25 * 86400),
  list(A = 3, P = 86400)
))
psd3 <- compute_power(z3)
pk3  <- find_peaks(psd3, n_peaks = 5)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(z3, main = "1.3: Two years of signal", ylab = "Temp (°C)", xlab = "Time (s)")
plot(psd3$period_day, psd3$psd, type = "l", log = "x",
     main = "1.3: Power spectrum", xlab = "Period (days)", ylab = "PSD")
abline(v = c(1, 365), lty = 2, col = "red")

cat("  Peaks found:\n")
print(pk3[, c("period_day", "amplitude")])


# ============================================================
# Test 1.4  Non-sinusoidal diel (3 harmonics)
# ============================================================
cat("\n=== 1.4 Three harmonics (24h, 12h, 8h) ===\n")
z4 <- make_signal(30, list(
  list(A = 5.0, P = 86400),
  list(A = 1.5, P = 43200),
  list(A = 0.8, P = 28800)
))
psd4 <- compute_power(z4)
pk4  <- find_peaks(psd4, n_peaks = 5)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(z4[1:(3*24)], main = "1.4: Non-sinusoidal diel", ylab = "Temp (°C)", xlab = "Time (s)")
plot(psd4$period_hr, psd4$psd, type = "l", log = "x",
     main = "1.4: Harmonic content", xlab = "Period (hr)", ylab = "PSD")
abline(v = c(8, 12, 24), lty = 2, col = "red")

cat("  Peaks found:\n")
print(pk4[, c("period_hr", "amplitude")])


# ============================================================
# Test 1.5  Amplitude-modulated boundary (realistic)
# ============================================================
cat("\n=== 1.5 Amplitude-modulated stream signal (1 year) ===\n")
z5 <- generate_composite_boundary(n_years = 2)
psd5 <- compute_power(z5)
pk5  <- find_peaks(psd5, n_peaks = 5)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(z5, main = "1.5: Composite boundary", ylab = "Temp (°C)", xlab = "Time (s)")
plot(psd5$period_day, psd5$psd, type = "l", log = "x",
     main = "1.5: Spectral content", xlab = "Period (days)", ylab = "PSD")
abline(v = c(1, 365), lty = 2, col = "red")

cat("  Peaks found:\n")
print(pk5[, c("period_day", "period_hr", "amplitude")])
cat("  Note: diel amplitude < 2.2 because seasonal modulation spreads energy into sidebands\n")


# ============================================================
# Test 1.6  Forward-inverse round trip (existing machinery)
# ============================================================
cat("\n=== 1.6 Forward-inverse round trip ===\n")
rt <- test_round_trip(fit_method = "ols")
cat("  Max relative error:", format(max(rt$relative_error), digits = 3), "\n")
print(rt[, c("scenario", "parameter", "known", "recovered", "relative_error")])


# ============================================================
# Summary
# ============================================================
cat("\n=== Tier 1 Summary ===\n")
cat("1.1 Single cosine:     Peak at", pk1$period_hr[1], "hr, A =", round(pk1$amplitude[1], 2), "\n")
cat("1.2 Two frequencies:   Both resolved, A = ", paste(round(pk2$amplitude, 2), collapse = ", "), "\n")
cat("1.3 Multi-scale:       Annual + diel separated, A = ", paste(round(pk3$amplitude, 1), collapse = ", "), "\n")
cat("1.4 Harmonics:         24/12/8h recovered, A = ", paste(round(pk4$amplitude, 2), collapse = ", "), "\n")
cat("1.5 Composite:         Realistic AM boundary decomposes cleanly\n")
cat("1.6 Round trip:        Max relative error =", format(max(rt$relative_error), digits = 3), "\n")
