
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

cat("\n=== Script complete. Explore objects: fwd, inv_ols, inv_nls, comparison ===\n")
