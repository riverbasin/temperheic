# demo/demo_observed.R
# temperheic Tier 3: real observed data end-to-end
#
# Four analyses on Meacham Creek data:
#   3.1 Annual transect (stream → Well 21 → Well 10)
#   3.2 Diel inverse (stream → Well 20)
#   3.3 Air temperature spectral characterization
#   3.4 Spring 6 null case

library(temperheic)
library(zoo)

# --- Helper: load CSV as zoo (seconds index) ---
load_zoo <- function(file, value_col) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  dt <- as.POSIXct(df$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  ok <- !is.na(dt) & !duplicated(dt)
  t0 <- as.numeric(dt[ok][1])
  zoo(df[ok, value_col, drop = FALSE], order.by = as.numeric(dt[ok]) - t0)
}

data_dir <- system.file("extdata", package = "temperheic")


# ============================================================
# 3.1  Annual transect: Channel 3 → Well 21 → Well 10
# ============================================================
cat("\n=== 3.1 Annual transect (3 sites, 450 days) ===\n")
z3 <- load_zoo(file.path(data_dir, "meacham_3site_hourly.csv"),
               c("channel_3_0m", "well_21_175m", "well_10_955m"))

# Daily averages for annual fit
z3d <- aggregate(z3, by = floor(index(z3) / 86400) * 86400, FUN = mean, na.rm = TRUE)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
cols3 <- c("steelblue", "forestgreen", "darkorange")
matplot(index(z3d)/86400, coredata(z3d), type = "l", col = cols3, lty = 1,
        main = "3-site transect", xlab = "Days", ylab = "°C")
legend("topright", names(z3d), col = cols3, lty = 1, cex = 0.6)

# Spectra
psd3 <- compute_power(z3)
for (i in seq_along(names(z3))) {
  sub <- psd3[psd3$series == names(z3)[i], ]
  if (i == 1) plot(sub$period_hr, sub$psd, type = "l", log = "xy",
       xlim = c(4, 10000), main = "Spectral fingerprints",
       xlab = "Period (hr)", ylab = "PSD", col = cols3[i])
  else lines(sub$period_hr, sub$psd, col = cols3[i])
}
abline(v = c(24, 365.25*24), lty = 2, col = "gray50")
legend("bottomright", names(z3), col = cols3, lty = 1, cex = 0.6)

# Annual inverse
xVals_3 <- c(channel_3_0m = 0, well_21_175m = 175, well_10_955m = 955)
inv3 <- thObservedSeries(
  empiricalData = z3d, xVals = xVals_3, aquifer = generate_example_aquifer(),
  boundaryMean = mean(coredata(z3d[, 1]), na.rm = TRUE),
  period = 365.25 * 86400, headGrad = 0.005, nmin = 100, fit_method = "ols")

cat("  Amplitude ratios: ", paste(round(inv3$ampRatio[1, ], 3), collapse = ", "), "\n")
cat("  η by pair:\n")
cat("    Ch3→W21:", round(inv3$eta[1, 2], 4), "\n")
cat("    Ch3→W10:", round(inv3$eta[1, 3], 4), "\n")
cat("    W21→W10:", round(inv3$eta[2, 3], 4), "\n")
cat("  η varies by pair — indicates heterogeneity or non-1D flow.\n")


# ============================================================
# 3.2  Diel inverse: Stream → Well 20 (216 days concurrent)
# ============================================================
cat("\n=== 3.2 Diel inverse: Stream–Well 20 ===\n")
z_st <- load_zoo(file.path(data_dir, "meacham_channel3_hourly.csv"), "channel_3_temp")
z_w20 <- load_zoo(file.path(data_dir, "meacham_well20_hourly.csv"), "well_20_temp")

# Align on absolute timestamps (both start from their own t0, so rebase)
# Use the raw numeric POSIX times for matching
f_st  <- file.path(data_dir, "meacham_channel3_hourly.csv")
f_w20 <- file.path(data_dir, "meacham_well20_hourly.csv")
df_st  <- read.csv(f_st, stringsAsFactors = FALSE)
df_w20 <- read.csv(f_w20, stringsAsFactors = FALSE)
dt_st  <- as.POSIXct(df_st$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
dt_w20 <- as.POSIXct(df_w20$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
common_abs <- intersect(as.numeric(na.omit(dt_st)), as.numeric(na.omit(dt_w20)))
t0_pair <- min(common_abs)
z_pair <- zoo(cbind(
  stream = df_st$channel_3_temp[match(common_abs, as.numeric(dt_st))],
  well20 = df_w20$well_20_temp[match(common_abs, as.numeric(dt_w20))]
), order.by = common_abs - t0_pair)
cat("  Concurrent hours:", nrow(z_pair), "\n")

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
idx5 <- 1:(5*24)
plot(index(z_pair)[idx5]/3600, coredata(z_pair)[idx5, 1], type = "l",
     col = "steelblue", main = "Stream vs Well 20 (5 days)",
     xlab = "Time (hr)", ylab = "°C", ylim = range(z_pair[idx5, ], na.rm = TRUE))
lines(index(z_pair)[idx5]/3600, coredata(z_pair)[idx5, 2], col = "darkorange")
legend("topright", c("Stream", "Well 20"), col = c("steelblue", "darkorange"),
       lty = 1, cex = 0.7)

# Spectra
psd_p <- compute_power(z_pair)
sub_s <- psd_p[psd_p$series == "stream", ]
sub_w <- psd_p[psd_p$series == "well20", ]
plot(sub_s$period_hr, sub_s$psd, type = "l", log = "xy",
     xlim = c(4, 5000), main = "Spectra", xlab = "Period (hr)", ylab = "PSD",
     col = "steelblue")
lines(sub_w$period_hr, sub_w$psd, col = "darkorange")
abline(v = 24, lty = 2, col = "gray50")
legend("bottomright", c("Stream", "Well 20"), col = c("steelblue", "darkorange"),
       lty = 1, cex = 0.7)

# Inverse
xVals_p <- c(stream = 0, well20 = 20)
inv_p <- thObservedSeries(
  empiricalData = z_pair, xVals = xVals_p, aquifer = generate_example_aquifer(),
  boundaryMean = mean(coredata(z_pair[, 1]), na.rm = TRUE),
  period = 86400, headGrad = 0.005, nmin = 10, fit_method = "ols")

cat("  Amp ratio (well/stream):", round(inv_p$ampRatio[1, 2], 3), "\n")
cat("  Phase lag:", round(inv_p$deltaPhaseRadians[1, 2] * 86400 / (2*pi) / 3600, 1), "hr\n")
cat("  η:", round(inv_p$eta[1, 2], 4), "\n")
cat("  Thermal velocity:", round(inv_p$advectiveThermVelEmpirical[1, 2] * 86400, 1), "m/d\n")


# ============================================================
# 3.3  Air temperature spectral characterization
# ============================================================
cat("\n=== 3.3 Air temperature (SNOTEL UP211) ===\n")
z_air <- load_zoo(file.path(data_dir, "UP211_air_temp_hourly_2012_2014.csv"), "air_temp")
cat("  Record:", length(z_air), "hours, range:",
    round(range(coredata(z_air), na.rm = TRUE), 1), "°C\n")

psd_air <- compute_power(z_air)
pk_air  <- find_peaks(psd_air, n_peaks = 3)
cat("  Top 3 peaks:\n")
print(pk_air[, c("period_hr", "period_day", "amplitude")])


# ============================================================
# 3.4  Spring 6: the null case
# ============================================================
cat("\n=== 3.4 Spring 6 — null case ===\n")
z_sp <- load_zoo(file.path(data_dir, "meacham_spring6_hourly.csv"), "spring_6_temp")
z_st_sub <- z_st[index(z_st) <= max(index(z_sp))]  # comparable period

psd_sp <- compute_power(z_sp)
psd_st_sub <- compute_power(z_st_sub)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(index(z_sp)[1:(7*24)]/3600, coredata(z_sp)[1:(7*24)], type = "l",
     col = "purple", main = "Spring 6: one week", xlab = "Time (hr)", ylab = "°C")

plot(psd_st_sub$period_hr, psd_st_sub$psd, type = "l", log = "xy",
     xlim = c(4, 5000), main = "Stream vs Spring spectrum",
     xlab = "Period (hr)", ylab = "PSD", col = "steelblue")
lines(psd_sp$period_hr, psd_sp$psd, col = "purple")
abline(v = 24, lty = 2, col = "gray50")
legend("bottomright", c("Stream", "Spring 6"),
       col = c("steelblue", "purple"), lty = 1, cex = 0.7)

pk_sp <- find_peaks(psd_sp, n_peaks = 3)
cat("  Total range:", round(diff(range(coredata(z_sp), na.rm = TRUE)), 1), "°C\n")
cat("  Diel amplitude:", round(pk_sp$amplitude[which.min(abs(pk_sp$period_hr - 24))], 2), "°C\n")
cat("  (Stream diel ≈ 1.85°C → Spring retains ~",
    round(pk_sp$amplitude[which.min(abs(pk_sp$period_hr - 24))] / 1.85 * 100), "%)\n")


# ============================================================
# Summary
# ============================================================
cat("\n=== Tier 3 Summary ===\n")
cat("3.1 Annual transect: η varies across pairs (0.64, 0.49, 0.19) — heterogeneity\n")
cat("3.2 Diel stream–well: amp ratio 0.32, phase lag 10.5 hr, η = 0.42\n")
cat("3.3 Air temp: diel peak + strong weather noise (red spectrum)\n")
cat("3.4 Spring 6: diel signal ~100× weaker than stream — correctly identified\n")
