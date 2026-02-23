# thSeries.R - Forward and inverse models for temperheic
#
# Contains:
#   thSeries()         - Forward model: generate synthetic temperature time series
#   thObservedSeries() - Inverse model: recover flux/diffusivity from observed data


#' Generate synthetic temperature time series (forward model)
#'
#' Given a thSignal object (which encapsulates aquifer, hydrology, and boundary
#' conditions), generates temperature time series at specified depths using the
#' analytical solution to the 1D advection-diffusion equation (Luce et al. 2013,
#' Eq. 6):
#'
#'   T(z,t) = T_mean + A * exp(-z / thermDecayDist) *
#'            cos((2*pi/P) * (t - phase - z/phaseVel))
#'
#' Also computes pairwise amplitude ratios, phase differences, depth
#' differences, and eta values for all sensor combinations.
#'
#' @param signal A thSignal object containing transport and boundary parameters.
#' @param xVals Numeric vector of sensor depths (m). Must include 0 (boundary).
#' @param tVals Numeric vector of time values (s) at which to evaluate.
#' @param specificUnits A thUnits object. Must match the signal units.
#'
#' @return A temperheic S3 object (class thSpecifiedSeries/thSeries/temperheic)
#'
#' @references
#' Stallman, R. W. (1965). Steady one-dimensional fluid flow in a semi-infinite
#' porous medium with sinusoidal surface temperature. Journal of Geophysical
#' Research, 70(12), 2821-2827.
#'
#' Luce, C. H., Tonina, D., Gariglio, F., & Applebee, R. (2013). Solutions for
#' the diurnally forced advection-diffusion equation to estimate bulk fluid
#' velocity and diffusivity in streambeds from temperature time series. Water
#' Resources Research, 49, 488-506. doi:10.1002/wrcr.20090
#' @export
thSeries <- function(signal, xVals, tVals, specificUnits) {

  if (!identical(attr(signal, "specificUnits"), specificUnits)) {
    stop("The units of the signal are not the same as those specified in specificUnits.")
  }

  sensor_names <- paste0("x", xVals)

  # --- Generate temperature at each depth via analytical solution ---
  # Luce 2013 Eq. 6:
  #   T(z,t) = mean + A*exp(-z/thermDecayDist)*cos(omega*(t - phase) - z*omega/phaseVel)
  timeSeries <- zoo::zoo(
    do.call(data.frame, args = purrr::map(xVals, function(x) {
      signal$boundary$mean +
        signal$boundary$amplitude *
        exp(-x / signal$thermDecayDist) *
        cos((2 * pi / signal$boundary$period) *
              ((tVals - tVals[1]) - signal$boundary$phase - (x / signal$phaseVel)))
    })),
    order.by = tVals
  )
  names(timeSeries) <- sensor_names

  # --- Amplitude at each depth: A(z) = A_surface * exp(-z/thermDecayDist) ---
  amplitude <- purrr::map_dbl(xVals, function(x) {
    signal$boundary$amplitude * exp(-x / signal$thermDecayDist)
  })
  names(amplitude) <- sensor_names

  # --- Phase at each depth (radians) ---
  phaseRadians <- purrr::map_dbl(xVals, function(x) {
    (x / signal$phaseVel + signal$boundary$phase) * (2 * pi / signal$boundary$period)
  })
  names(phaseRadians) <- sensor_names

  # --- Pairwise matrices for all sensor combinations ---
  combos    <- expand.grid(from = sensor_names, to = sensor_names)
  xVals_named <- stats::setNames(xVals, sensor_names)

  deltaXvals        <- derived2DArray(xVals_named[combos$to] - xVals_named[combos$from], sensor_names)
  ampRatio          <- derived2DArray(amplitude[combos$to] / amplitude[combos$from], sensor_names)
  deltaPhaseRadians <- derived2DArray(phaseRadians[combos$to] - phaseRadians[combos$from], sensor_names)
  eta               <- -log(ampRatio) / deltaPhaseRadians  # Luce et al. 2013

  structure(
    list(
      # User inputs
      signal         = signal,
      xVals          = xVals_named,
      # Derived: per-sensor
      timeSeries     = timeSeries,
      amplitude      = amplitude,
      phaseRadians   = phaseRadians,
      # Derived: pairwise matrices
      ampRatio          = ampRatio,
      deltaPhaseRadians = deltaPhaseRadians,
      deltaXvals        = deltaXvals,
      eta               = eta
    ),
    class = c("thSpecifiedSeries", "thSeries", "temperheic"),
    units = .thSpecificUnits("T", specificUnits),
    derivedValueNames = c("timeSeries", "amplitude", "phaseRadians",
                          "ampRatio", "deltaPhaseRadians", "deltaXvals", "eta"),
    specificUnits = specificUnits,
    thObjectNames = "signal"
  )
}


#' Recover thermal parameters from observed temperature data (inverse model)
#'
#' Fits cosine curves to multi-depth temperature time series and uses the
#' amplitude ratios and phase differences to estimate thermal velocity,
#' diffusivity, and Darcy flux via the eta-parameter method (Luce et al. 2013).
#'
#' The estimation pipeline:
#' 1. fit_ols() or fitCosine() extracts amplitude and phase at each sensor depth
#' 2. Pairwise amplitude ratios and phase differences are computed
#' 3. eta = -ln(Ar) / delta_phi  (Luce 2013)
#' 4. Thermal velocity from eta and phase diff (Luce et al. 2013, thermal velocity from eta)
#' 5. Effective diffusivity from eta and amplitude ratio (Luce et al. 2013, diffusivity from eta)
#' 6. Darcy flux, water velocity, dispersivity, hydraulic conductivity derived
#'
#' @param empiricalData A zoo object with one column per sensor depth.
#' @param xVals Named numeric vector of sensor depths (m). Must include 0.
#'   Names must match column names of empiricalData.
#' @param aquifer A thAquifer object with thermal/physical properties.
#' @param boundaryMean Numeric - mean temperature at the boundary (degC).
#' @param period Numeric - period of the target cycle (s). 86400 for diel.
#' @param headGrad Numeric - hydraulic head gradient (dimensionless).
#' @param nmin Integer - minimum non-NA observations to attempt fitting.
#' @param freq Numeric - angular frequency (rad/s). Default 2*pi/period.
#' @param optimizeRange Numeric length-2 - phase window as fractions of period.
#' @param specificUnits A thUnits object. Must match the aquifer units.
#' @param empiricalDataPeriods Numeric vector - number of complete cycles per
#'   sensor. Default rep(1, ncol(empiricalData)).
#' @param fit_method Character -- method for amplitude/phase extraction.
#'   "ols" (default) uses linearized OLS harmonic regression (Luce et al. 2013).
#'   "nls" uses the legacy nonlinear least squares cosine fit (fitCosine).
#'   "fft" uses FFT spectral extraction at the nearest Fourier bin.
#'
#' @return A temperheic S3 object (class thObservedSeries/thSeries/temperheic)
#'
#' @references
#' Luce, C. H., Tonina, D., Gariglio, F., & Applebee, R. (2013). Solutions for
#' the diurnally forced advection-diffusion equation to estimate bulk fluid
#' velocity and diffusivity in streambeds from temperature time series. Water
#' Resources Research, 49, 488-506. doi:10.1002/wrcr.20090
#'
#' Luce, C. H., & Tonina, D. (2017). Scaling thermal transport from streams
#' and their beds. Water Resources Research, 53, 9771-9790.
#' doi:10.1002/2017WR021472
#' @export
thObservedSeries <- function(empiricalData,
                             xVals,
                             aquifer,
                             boundaryMean,
                             period,
                             headGrad,
                             nmin,
                             freq = (2 * pi) / period,
                             optimizeRange = c(-1/8, 7/8),
                             specificUnits = thUnits(),
                             empiricalDataPeriods = rep(1, ncol(empiricalData)),
                             fit_method = c("ols", "nls", "fft")) {

  # --- Input validation ---
  fit_method <- match.arg(fit_method)

  if ((optimizeRange[2] - optimizeRange[1]) != 1) {
    stop("max optimize range - min optimize range must = 1.0")
  }
  if (!inherits(empiricalData, "zoo")) {
    stop("empiricalData must be a zoo object.")
  }
  if (!identical(sort(names(empiricalData)), sort(names(xVals)))) {
    stop("The names of the columns in empiricalData must match the names of xVals.")
  }
  if (length(xVals) != ncol(empiricalData)) {
    stop("Please specify one distance for each column in empiricalData.")
  }
  if (!identical(specificUnits, attr(aquifer, "specificUnits"))) {
    stop("specificUnits of aquifer does not match the specificUnits argument.")
  }

  # Ensure ascending depth order; boundary (depth=0) must be present
  xVals <- xVals[order(xVals)]
  if (xVals[1] != 0) stop("xVals must include zero to designate input signal.")
  empiricalData <- empiricalData[, names(xVals)]
  sensor_names <- names(xVals)
  nSeries <- ncol(empiricalData)

  # --- Fit cosine to each sensor and get pairwise comparisons ---
  # Dispatch to the selected fitting method. All methods return the same
  # output structure: list(deltaPhaseRadians, ampRatio) with attributes
  # "amplitudes" and "phases".
  fit_fn <- switch(fit_method,
    ols = fit_ols,
    nls = fitCosine,
    fft = fit_fft,
    stop("Unknown fit_method: ", fit_method)
  )
  results <- fit_fn(empiricalData, boundaryMean, period,
                    optimizeRange, nmin, empiricalDataPeriods)

  amplitude     <- attr(results, "amplitudes")
  relativePhase <- attr(results, "phases")

  # --- Build pairwise matrices ---
  combos <- expand.grid(from = sensor_names, to = sensor_names)
  diag_locs <- (0:(nSeries - 1)) * nSeries + 1:nSeries

  ampRatio          <- derived2DArray(results$ampRatio, sensor_names)
  deltaPhaseRadians <- derived2DArray(results$deltaPhaseRadians, sensor_names)
  deltaXvals        <- derived2DArray(xVals[combos$to] - xVals[combos$from], sensor_names)

  # eta = -ln(Ar) / delta_phi  (Luce et al. 2013)
  eta_vec <- -log(results$ampRatio) / results$deltaPhaseRadians
  eta_vec[diag_locs] <- NaN  # self-comparisons undefined
  eta <- derived2DArray(eta_vec, sensor_names)

  # --- Thermal velocity (Luce et al. 2013; v_t from eta and delta_phi) ---
  # v_t = (omega * dz / sqrt(ln(Ar)^2 + dphi^2)) * (1 - eta^2) / sqrt(1 + eta^2)
  advectiveThermVel <- ((freq * deltaXvals) /
    (sqrt(log(ampRatio)^2 + deltaPhaseRadians^2))) *
    ((1 - eta^2) / sqrt(1 + eta^2))

  # --- Effective diffusivity (Luce et al. 2013; kappa_e from eta and Ar) ---
  # kappa_e = eta * omega * dz^2 / (ln(Ar)^2 + dphi^2)
  diffusivity_effective <- (eta * freq * deltaXvals^2) /
    (log(ampRatio)^2 + deltaPhaseRadians^2)

  # --- Derived physical parameters ---
  # Darcy flux: q = v_t * (rho_m*c_m) / (rho_w*c_w)
  darcyFlux <- advectiveThermVel *
    (aquifer$volHeatCap_bulk / (aquifer$density_h2o * aquifer$spHeat_h2o))

  velocity_h2o <- darcyFlux / aquifer$porosity

  # Dispersivity: psi = (D*rho_m*c_m - lambda_m) / (v_t * rho_m*c_m)
  dispersivity <- ((diffusivity_effective * aquifer$volHeatCap_bulk -
                      aquifer$thermCond_bulk) /
                     (advectiveThermVel * aquifer$volHeatCap_bulk))

  hydraulicCond <- darcyFlux / headGrad

  # --- Thermal decay distance and phase velocity ---
  # These use the full expressions from the PDE solution (Luce 2013 Eqs. 8-9)
  rad1  <- (advectiveThermVel^4 + (4 * 2 * pi * freq * diffusivity_effective)^2)^0.5

  thermDecayDist <- 2 * diffusivity_effective *
    (sqrt(2) / (sqrt(rad1 + advectiveThermVel^2) - sqrt(2) * advectiveThermVel))

  phaseVel <- 2 * diffusivity_effective *
    (sqrt(2) / abs(sqrt(rad1 - advectiveThermVel^2))) * 2 * pi * freq

  pecletNumber <- (advectiveThermVel * thermDecayDist) / diffusivity_effective

  # --- Build explicit object ---
  # Only include meaningful elements; exclude config params and intermediates
  structure(
    list(
      # User inputs
      empiricalData = empiricalData,
      xVals         = xVals,
      aquifer       = aquifer,
      boundaryMean  = boundaryMean,
      period        = period,
      headGrad      = headGrad,
      freq          = freq,
      # Derived: per-sensor
      amplitude     = amplitude,
      relativePhase = relativePhase,
      # Derived: pairwise matrices (core eta-method results)
      ampRatio          = ampRatio,
      deltaPhaseRadians = deltaPhaseRadians,
      eta               = eta,
      # Derived: physical parameter estimates (pairwise matrices)
      advectiveThermVelEmpirical      = advectiveThermVel,
      diffusivity_effective_empirical = diffusivity_effective,
      darcyFlux                       = darcyFlux,
      velocity_h2o                    = velocity_h2o,
      dispersivity                    = dispersivity,
      hydraulicCond                   = hydraulicCond,
      thermDecayDist                  = thermDecayDist,
      phaseVel                        = phaseVel,
      pecletNumber                    = pecletNumber
    ),
    class = c("thObservedSeries", "thSeries", "temperheic"),
    units = .thSpecificUnits("T", specificUnits),
    derivedValueNames = c("amplitude", "relativePhase",
                          "ampRatio", "deltaPhaseRadians", "eta",
                          "advectiveThermVelEmpirical",
                          "diffusivity_effective_empirical",
                          "darcyFlux", "velocity_h2o", "dispersivity",
                          "hydraulicCond", "thermDecayDist", "phaseVel",
                          "pecletNumber"),
    specificUnits = specificUnits,
    thObjectNames = "aquifer"
  )
}
