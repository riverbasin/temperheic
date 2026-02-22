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
#'
#' @return A list mirroring the structure of series_list, with thObservedSeries
#'   objects replacing thSeries objects.
#' @export
generate_example_observed <- function(
    series_list,
    nmin = list(diel = 10, annual = 180)) {

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
        specificUnits = thUnits()
      )
    })
  })
}


#' Round-trip validation: forward -> inverse -> compare
#'
#' Generates synthetic data with known parameters, runs the inverse model,
#' and compares recovered values to ground truth. Returns a tidy data frame
#' with relative errors for key parameters.
#'
#' @param aquifer A thAquifer object. Default uses generate_example_aquifer().
#' @param scenarios Passed to generate_example_series(). Default scenarios.
#'
#' @return A tibble with columns: scenario, combo, parameter, known, recovered,
#'   relative_error. Uses the [1,2] element of pairwise matrices (boundary
#'   sensor vs. first subsurface sensor).
#' @export
test_round_trip <- function(
    aquifer = generate_example_aquifer(),
    scenarios = NULL) {

  # Generate forward and inverse
  if (is.null(scenarios)) {
    fwd_list <- generate_example_series(aquifer)
  } else {
    fwd_list <- generate_example_series(aquifer, scenarios)
  }
  inv_list <- generate_example_observed(fwd_list)

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
