# th_test.R -- Round-trip validation and example generation for temperheic
#
# Provides helpers to:
#   1. Create example objects with realistic parameter sets
#   2. Run forward -> inverse round-trip validation
#   3. Exercise parameter combinations for sensitivity testing
#
# These are development/testing utilities, not end-user analysis functions.


#' Create an example thAquifer with basalt-like defaults
#'
#' Default values from Columbia River basalt (Burns et al. 2015, GeoFluids v.15).
#'
#' @return A thAquifer object.
#'
#' @references
#' Burns, E. R., Williams, C. F., Ingebritsen, S. E., Voss, C. I., Spane,
#' F. A., & DeAngelo, J. (2015). Understanding heat and groundwater flow
#' through continental flood basalt provinces: insights gained from alternative
#' models of permeability/depth relationships for the Columbia Plateau, USA.
#' Geofluids, 15, 120-138.
#' @export
generate_example_aquifer <- function(
    porosity = 0.25,
    thermCond_sed = 0.0016,
    thermCond_h2o = 0.000598,
    spHeat_sed = 0.84,
    spHeat_h2o = 4.186,
    density_sed = 3000,
    density_h2o = 998) {
  thAquifer(porosity, thermCond_sed, thermCond_h2o,
            spHeat_sed, spHeat_h2o, density_sed, density_h2o)
}


#' Generate forward-model time series across parameter combinations
#'
#' Creates thSeries objects for factorial combinations of hydraulic
#' conductivity and dispersivity, at both diel and annual time scales.
#'
#' @param aquifer A thAquifer object. Default uses generate_example_aquifer().
#' @param scenarios A named list of scenario specifications.
#'
#' @return A named list of lists. Top level is scenario name (e.g., 'diel',
#'   'annual'). Each contains one thSeries per parameter combination.
#' @export
generate_example_series <- function(
    aquifer = generate_example_aquifer(),
    scenarios = list(
      diel = list(
        k           = c(10) / 86400,
        dispersivity = c(0.001, 0.01, 0.1),
        headGrad    = 0.1,
        mean        = 22,
        amplitude   = 6,
        phase       = 86400 + 3600,
        period      = 86400,
        xVals       = c(0, 0.33, 0.66, 1),
        tVals       = seq(0, 3600 * 48, 3600)
      ),
      annual = list(
        k           = c(100) / 86400,
        dispersivity = c(100, 1000),
        headGrad    = 0.01,
        mean        = 11.78,
        amplitude   = 8.87,
        phase       = 1702226,
        period      = 86400 * 365,
        xVals       = c(0, 120, 240, 360),
        tVals       = seq(0, 2 * 365 * 86400, 86400)
      )
    )) {

  purrr::map(scenarios, function(sc) {
    # Build parameter grid: all combos of k and dispersivity
    params <- expand.grid(
      k    = sc$k,
      disp = sc$dispersivity,
      stringsAsFactors = FALSE
    )

    # One thSeries per parameter combination
    purrr::pmap(params, function(k, disp) {
      hy  <- thHydro(k, disp, sc$headGrad, aquifer)
      bo  <- thBoundary(sc$mean, sc$amplitude, sc$phase, sc$period)
      sig <- thSignal(hy, bo)
      thSeries(sig, xVals = sc$xVals, tVals = sc$tVals, specificUnits = thUnits())
    }) |> stats::setNames(paste0("k", params$k, "_disp", params$disp))
  })
}


#' Run inverse model on forward-model output
#'
#' Takes the output of generate_example_series() and runs thObservedSeries
#' on each, recovering thermal parameters from the synthetic temperature data.
#'
#' @param series_list Output from generate_example_series().
#' @param nmin Named list of minimum non-NA obs per sensor by scenario.
#' @param fit_method Character -- fitting method passed to thObservedSeries.
#'
#' @return A list mirroring the structure of series_list, with thObservedSeries
#'   objects replacing thSeries objects.
#' @export
generate_example_observed <- function(
    series_list,
    nmin = list(diel = 10, annual = 180),
    fit_method = "ols") {

  purrr::imap(series_list, function(scenario_series, scenario_name) {
    n <- nmin[[scenario_name]] %||% 10
    purrr::map(scenario_series, function(ser) {
      thObservedSeries(
        empiricalData = ser$timeSeries,
        xVals         = ser$xVals,
        aquifer       = ser$signal$hydro$aquifer,
        boundaryMean  = ser$signal$boundary$mean,
        period        = ser$signal$boundary$period,
        headGrad      = ser$signal$hydro$headGrad,
        nmin          = n,
        specificUnits = thUnits(),
        fit_method    = fit_method
      )
    })
  })
}


# Suppress R CMD check NOTEs for tibble column names used in NSE context
utils::globalVariables(c("known", "recovered"))

#' Round-trip validation: forward -> inverse -> compare
#'
#' Generates synthetic data with known parameters, runs the inverse model,
#' and compares recovered values to ground truth. Returns a tidy data frame
#' with relative errors for key parameters.
#'
#' @param aquifer A thAquifer object. Default uses generate_example_aquifer().
#' @param scenarios Passed to generate_example_series(). Default scenarios.
#' @param fit_method Character -- fitting method: "ols" (default), "nls", or "fft".
#'
#' @return A tibble with columns: scenario, combo, parameter, known, recovered,
#'   relative_error. Uses the [1,2] element of pairwise matrices (boundary
#'   sensor vs. first subsurface sensor).
#' @export
test_round_trip <- function(
    aquifer = generate_example_aquifer(),
    scenarios = NULL,
    fit_method = "ols") {

  # Generate forward and inverse
  if (is.null(scenarios)) {
    fwd_list <- generate_example_series(aquifer)
  } else {
    fwd_list <- generate_example_series(aquifer, scenarios)
  }
  inv_list <- generate_example_observed(fwd_list, fit_method = fit_method)

  # Compare known vs recovered for each combination
  purrr::imap_dfr(fwd_list, function(scenario_fwd, scenario_name) {
    scenario_inv <- inv_list[[scenario_name]]

    purrr::imap_dfr(scenario_fwd, function(fwd, combo_name) {
      inv <- scenario_inv[[combo_name]]

      # Known values from the signal object
      known_darcy <- fwd$signal$hydro$darcyFlux
      known_diff  <- fwd$signal$hydro$diffusivity_effective
      known_vt    <- fwd$signal$hydro$advectiveThermVel

      # Recovered values: [1,2] = boundary vs first subsurface sensor
      rec_darcy <- inv$darcyFlux[1, 2]
      rec_diff  <- inv$diffusivity_effective_empirical[1, 2]
      rec_vt    <- inv$advectiveThermVelEmpirical[1, 2]

      tibble::tibble(
        scenario = scenario_name,
        combo    = combo_name,
        parameter = c("darcyFlux", "diffusivity_effective", "advectiveThermVel"),
        known     = c(known_darcy, known_diff, known_vt),
        recovered = c(rec_darcy, rec_diff, rec_vt),
        relative_error = abs(recovered - known) / abs(known)
      )
    })
  })
}


#' Generate a composite stream temperature signal with diel and annual components
#'
#' Creates a synthetic hourly temperature time series that mimics real stream
#' temperature data, with amplitude-modulated diel cycles nested within a
#' seasonal (annual) envelope. Parameters are derived from Meacham Creek
#' Channel 3 observations.
#'
#' The signal model is a diel cosine whose amplitude is seasonally modulated
#' (larger in summer, smaller in winter), superimposed on an annual cosine.
#' This creates a realistic amplitude-modulated signal whose FFT contains
#' the fundamental diel and annual peaks plus sidebands at 1/day +/- 1/year.
#'
#' @param n_years Numeric, duration in years (default 1).
#' @param dt Numeric, time step in seconds (default 3600 = hourly).
#' @param T_mu Mean annual temperature (default 10.5 degC).
#' @param A_annual Annual amplitude (default 8.0 degC).
#' @param A_diel_mean Mean diel amplitude (default 2.2 degC).
#' @param A_diel_mod Seasonal modulation of diel amplitude (default 1.8 degC).
#' @param phi_diel_hour Hour of daily peak temperature (default 15 = 3 PM).
#' @param phi_annual_day Day of year for annual peak (default 200, late July).
#'
#' @return A zoo object with time index in seconds from start, suitable
#'   for use with fit_ols, fit_fft, or as boundary condition input.
#'
#' @references
#' Parameters derived from Meacham Creek Channel 3 (2012-2014).
#' @export
generate_composite_boundary <- function(
    n_years = 1,
    dt = 3600,
    T_mu = 10.5,
    A_annual = 8.0,
    A_diel_mean = 2.2,
    A_diel_mod = 1.8,
    phi_diel_hour = 15,
    phi_annual_day = 200
) {

  n_seconds <- round(n_years * 365.25 * 86400)
  t_sec <- seq(0, n_seconds - dt, by = dt)

  omega_annual <- 2 * pi / (365.25 * 86400)
  omega_diel   <- 2 * pi / 86400

  phi_annual   <- -omega_annual * phi_annual_day * 86400
  phi_diel     <- -omega_diel * phi_diel_hour * 3600
  phi_diel_mod <- phi_annual

  annual_component <- A_annual * cos(omega_annual * t_sec + phi_annual)
  diel_envelope    <- A_diel_mean + A_diel_mod * cos(omega_annual * t_sec + phi_diel_mod)
  diel_component   <- diel_envelope * cos(omega_diel * t_sec + phi_diel)

  temperature <- T_mu + annual_component + diel_component

  zoo::zoo(temperature, order.by = t_sec)
}
