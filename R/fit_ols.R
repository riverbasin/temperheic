
#' Fit amplitude and phase using OLS harmonic regression
#'
#' Uses the trigonometric identity A*cos(wt + phi) = alpha_s*sin(wt) + alpha_c*cos(wt)
#' to linearize cosine fitting into a standard lm() call (Luce et al. 2013;
#' observedAmpPhase.tex Eqs. 1-3). This avoids the convergence issues,
#' starting-value sensitivity, and noise-injection workarounds required by nls.
#'
#' The model at each depth is:
#'   T(t) = beta_0 + alpha_s * sin(omega * t) + alpha_c * cos(omega * t)
#'
#' where beta_0 captures the mean, and amplitude and phase are recovered via:
#'   A   = sqrt(alpha_s^2 + alpha_c^2)
#'   phi = atan2(alpha_s, alpha_c)
#'
#' @details
#' **Why OLS over nls:** The linearized form reduces cosine fitting to ordinary
#' least squares, which has a closed-form solution. No starting values, no
#' iteration, no convergence failure on perfect cosines. On synthetic round-trip
#' tests this recovers parameters to machine precision (~1e-15 relative error)
#' versus ~1e-4 for the nls approach with its required noise injection.
#'
#' **Phase convention:** atan2(alpha_s, alpha_c) returns the phase in (-pi, pi].
#' We wrap to [0, 2*pi) for consistency, then unwrap relative to the boundary
#' sensor phase. This is equivalent to the conditional arctan in
#' observedAmpPhase.tex Eq. 3 but handles all four quadrants automatically.
#'
#' **Interface compatibility:** Output structure matches fitCosine() exactly,
#' so thObservedSeries needs no changes to use this as a drop-in replacement.
#'
#' @param empiricalData A zoo object with one column per sensor depth.
#' @param boundaryMean Numeric -- mean temperature. Not used by OLS (beta_0 is
#'   estimated freely) but retained for interface compatibility with fitCosine.
#' @param periodInSeconds Numeric -- period of the target cycle (s).
#'   86400 for diel, 86400*365.25 for annual.
#' @param optimizeRange Numeric length-2 -- valid phase window as fractions of
#'   period relative to boundary sensor. Default c(-0.05, 1.05).
#' @param nmin Integer -- minimum non-NA observations to attempt fitting.
#' @param empiricalDataPeriods Numeric -- number of complete cycles in the data.
#'
#' @return A list with components:
#'   \describe{
#'     \item{deltaPhaseRadians}{Numeric vector of length n^2 -- pairwise phase
#'       differences (radians). Reshape with derived2DArray().}
#'     \item{ampRatio}{Numeric vector of length n^2 -- pairwise amplitude ratios.
#'       Reshape with derived2DArray().}
#'   }
#'   Attributes:
#'   \describe{
#'     \item{amplitudes}{Named numeric vector -- amplitude per sensor.}
#'     \item{phases}{Named numeric vector -- phase per sensor (radians).}
#'   }
#'
#' @references
#' Luce, C. H., Tonina, D., Gariglio, F., & Applebee, R. (2013). Solutions for
#' the diurnally forced advection-diffusion equation to estimate bulk fluid
#' velocity and diffusivity in streambeds from temperature time series. Water
#' Resources Research, 49, 488-506.
fit_ols <- function(empiricalData, boundaryMean, periodInSeconds,
                    optimizeRange, nmin, empiricalDataPeriods) {

  omega <- 2 * pi / periodInSeconds

  # Time vector in seconds from start (lm needs numeric, not POSIXct)
  seconds <- as.numeric(zoo::index(empiricalData))
  obsTime <- seconds - seconds[1]

  # Precompute basis vectors once -- shared across all sensors
  sin_basis <- sin(omega * obsTime)
  cos_basis <- cos(omega * obsTime)

  # --- Fit OLS to each sensor independently ---
  #
  # Model: T(t) = beta_0 + alpha_s * sin(omega*t) + alpha_c * cos(omega*t)
  # lm() gives closed-form solution: no iteration, no starting values.
  fits <- purrr::map(as.data.frame(empiricalData), function(eD) {

    # Require minimum non-NA observations
    if (length(stats::na.omit(eD)) < nmin) {
      return(c(fitAmp = NA_real_, fitPhase = NA_real_))
    }

    fit <- stats::lm(eD ~ sin_basis + cos_basis, na.action = stats::na.exclude)

    alpha_s <- stats::coef(fit)[["sin_basis"]]
    alpha_c <- stats::coef(fit)[["cos_basis"]]

    # Recover amplitude: A = sqrt(alpha_s^2 + alpha_c^2)
    amp <- sqrt(alpha_s^2 + alpha_c^2)

    # Recover phase: atan2 handles all four quadrants automatically,
    # equivalent to the conditional arctan in observedAmpPhase.tex Eq. 3
    pha <- atan2(alpha_s, alpha_c)

    # Wrap to [0, 2*pi)
    pha <- pha %% (2 * pi)

    c(fitAmp = amp, fitPhase = pha)
  })

  # Extract per-sensor amplitude and phase vectors
  fitAmp   <- purrr::map_dbl(fits, "fitAmp")
  fitPhase <- purrr::map_dbl(fits, "fitPhase")

  # --- Phase unwrapping relative to boundary sensor ---
  #
  # Boundary (first sensor, shallowest depth) sets the reference phase.

  # Deeper sensors should fall within [initialPhase + optimizeRange[1]*2pi,
  # initialPhase + optimizeRange[2]*2pi]. If a phase falls below this window,
  # add 2*pi to bring it into range.
  initialPhase <- fitPhase[1]
  relativeRange <- 2 * pi * optimizeRange + initialPhase
  fitPhase[fitPhase < relativeRange[1]] <- fitPhase[fitPhase < relativeRange[1]] + 2 * pi
  fitPhase <- fitPhase + (empiricalDataPeriods - 1) * 2 * pi

  # --- Pairwise amplitude ratios and phase differences ---
  #
  # expand.grid produces all (from, to) index pairs. The resulting length-n^2
  # vectors are reshaped into n x n matrices by derived2DArray() in the caller.
  combos <- expand.grid(from = seq_along(fitPhase), to = seq_along(fitAmp))
  deltaPhaseRadians <- fitPhase[combos$to] - fitPhase[combos$from]
  ampRatio          <- fitAmp[combos$to]   / fitAmp[combos$from]

  structure(
    list(deltaPhaseRadians = deltaPhaseRadians, ampRatio = ampRatio),
    amplitudes = fitAmp,
    phases     = fitPhase
  )
}

