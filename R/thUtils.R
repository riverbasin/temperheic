# thUtils.R — Utility functions for temperheic
#
# Contains:
#   fitCosine()      — Cosine fitting engine for the inverse model (thObservedSeries)
#   derived2DArray() — Builds named pairwise matrices from vectors of sensor comparisons


#' Build a named square matrix from a pairwise comparison vector
#'
#' Takes a vector of length n^2 (from expand.grid permutations of n sensors)
#' and reshapes it into an n x n matrix with labeled rows and columns.
#' Used by thObservedSeries to organize amplitude ratios, phase differences,
#' eta values, and depth differences into sensor-pair matrices.
#'
#' @param x Numeric vector of length n^2 (pairwise values)
#' @param seriesNames Character vector of length n (sensor/depth names)
#' @return A named n x n matrix with dimnames "from" (rows) and "to" (columns)
derived2DArray <- function(x, seriesNames) {
  array(
    data = x,
    dim = c(length(seriesNames), length(seriesNames)),
    dimnames = list(from = seriesNames, to = seriesNames)
  )
}


#' Fit cosine curves to multi-depth temperature time series
#'
#' Uses nonlinear least squares (nls) to fit a cosine function to each sensor's
#' temperature data, extracting amplitude and phase. Returns pairwise amplitude
#' ratios and phase differences for all sensor combinations.
#'
#' The fitted model at each depth is:
#'   T(t) = boundaryMean + A * cos((2*pi/period) * t - phi)
#'
#' where A (amplitude) and phi (phase) are estimated by nls, and boundaryMean
#' and period are held fixed (known from the boundary condition).
#'
#' @details
#' **Negative amplitude handling:** nls may converge to a negative amplitude
#' with a phase offset. When this happens, the amplitude is negated and pi
#' radians is added to the phase (mathematically equivalent: -A*cos(t - phi) =
#' A*cos(t - phi - pi)).
#'
#' **Phase wrapping:** Phases are adjusted relative to the boundary sensor's
#' phase to ensure consistent phase differences across depths. The
#' `optimizeRange` parameter defines the valid window (in periods) for phase
#' values relative to the boundary, and `empiricalDataPeriods` accounts for
#' multi-period datasets where deeper sensors may lag by more than one cycle.
#'
#' **Noise injection for synthetic data:** When input data is a perfect cosine
#' (zero residuals), nls encounters a singular gradient and fails. A negligible
#' amount of noise (~0.01% of the mean) is added to prevent this. This has no
#' measurable effect on real data but enables round-trip validation with
#' synthetic signals.
#'
#' @param empiricalData A zoo object with one column per sensor depth.
#'   Column names should correspond to sensor identifiers.
#' @param boundaryMean Numeric scalar — mean temperature at the boundary (degC).
#'   Held fixed during fitting.
#' @param periodInSeconds Numeric scalar — period of the target cycle (s).
#'   86400 for diel, 86400*365.25 for annual.
#' @param optimizeRange Numeric vector of length 2 — valid phase window as
#'   fractions of the period relative to the boundary sensor phase.
#'   Default c(-0.05, 1.05) allows slight negative phases and up to one
#'   full period of lag.
#' @param nmin Integer — minimum number of non-NA observations required to
#'   attempt fitting for a given sensor. Sensors below this threshold return
#'   NA for amplitude and phase.
#' @param empiricalDataPeriods Numeric — number of complete cycles in the data.
#'   Used for phase unwrapping in multi-period datasets (typically 1 for
#'   standard analysis windows).
#'
#' @return A list with two components:
#'   \describe{
#'     \item{deltaPhaseRadians}{Numeric vector of length n^2 — pairwise phase
#'       differences (radians). Reshape with derived2DArray().}
#'     \item{ampRatio}{Numeric vector of length n^2 — pairwise amplitude ratios.
#'       Reshape with derived2DArray().}
#'   }
#'   Attributes:
#'   \describe{
#'     \item{amplitudes}{Named numeric vector — fitted amplitude per sensor.}
#'     \item{phases}{Named numeric vector — fitted phase per sensor (radians).}
#'   }
fitCosine <- function(empiricalData, boundaryMean, periodInSeconds,
                      optimizeRange, nmin, empiricalDataPeriods) {

  # Initial amplitude estimates (half the observed range at each depth)
  ampEst <- purrr::map_dbl(as.data.frame(empiricalData),
                            ~ (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)) / 2)

  # Time vector in seconds from start (nls needs numeric, not POSIXct)
  seconds <- as.numeric(zoo::index(empiricalData))
  obsTime <- seconds - seconds[1]

  # --- Fit cosine to each sensor independently via nls ---
  #
  # Each sensor gets: T(t) = boundaryMean + A * cos(omega*t - phi)
  # with A and phi as free parameters. boundaryMean and period are fixed.
  fits <- purrr::map2(
    as.data.frame(empiricalData),
    ampEst,
    function(eD, amp_init) {
      # Require minimum number of non-NA observations
      if (length(stats::na.omit(eD)) < nmin) {
        return(c(fitAmp = NA_real_, fitPhase = NA_real_))
      }

      # Inject negligible noise to avoid singular gradient on perfect cosines.
      # nls requires non-zero residuals to compute the Jacobian; synthetic data
      # from the forward model (thSeries) produces exact cosines that cause
      # convergence failure without this.
      noise_scale <- mean(eD, na.rm = TRUE) * 1e-4
      eD <- eD + stats::runif(length(eD), min = -noise_scale, max = noise_scale)

      fit <- stats::nls(
        formula = eD ~ boundaryMean + AmpY._ * cos((2 * pi / periodInSeconds) * obsTime - PhaY._),
        start   = list(AmpY._ = amp_init, PhaY._ = 0),
        na.action = "na.exclude"
      )

      amp <- stats::coef(fit)[["AmpY._"]]
      pha <- stats::coef(fit)[["PhaY._"]]

      # Normalize negative amplitudes: -A*cos(t - phi) == A*cos(t - phi - pi)
      if (amp < 0) {
        amp <- -amp
        pha <- pha + pi
      }

      # Wrap phase to [0, 2*pi)
      pha <- pha %% (2 * pi)

      c(fitAmp = amp, fitPhase = pha)
    }
  )

  # Extract amplitude and phase vectors from the list of 2-element results
  fits_mat <- do.call(rbind, fits)
  fitAmp   <- fits_mat[, "fitAmp"]
  fitPhase <- fits_mat[, "fitPhase"]
  names(fitAmp)   <- names(empiricalData)
  names(fitPhase) <- names(empiricalData)

  # --- Phase unwrapping relative to boundary sensor ---
  #
  # The boundary sensor (first column, shallowest depth) sets the reference
  # phase. Deeper sensors should have phases that fall within
  # [initialPhase + optimizeRange[1]*2*pi, initialPhase + optimizeRange[2]*2*pi].
  #
  # If a fitted phase falls below this window (e.g., nls found a phase near 0

  # when the boundary is near 2*pi), add 2*pi to bring it into range.
  #
  # empiricalDataPeriods accounts for datasets spanning multiple cycles --
  # deeper sensors may lag by more than one full period.
  initialPhase <- fitPhase[1]
  relativeRange <- 2 * pi * optimizeRange + initialPhase
  fitPhase[fitPhase < relativeRange[1]] <- fitPhase[fitPhase < relativeRange[1]] + 2 * pi
  fitPhase <- fitPhase + (empiricalDataPeriods - 1) * 2 * pi

  # --- Compute all pairwise amplitude ratios and phase differences ---
  #
  # expand.grid produces all (from, to) sensor index pairs. The resulting
  # vectors are length n^2 and can be reshaped into n x n matrices by
  # derived2DArray() in the calling function (thObservedSeries).
  combos <- expand.grid(from = seq_along(fitPhase), to = seq_along(fitAmp))
  deltaPhaseRadians <- fitPhase[combos$to] - fitPhase[combos$from]
  ampRatio          <- fitAmp[combos$to]   / fitAmp[combos$from]

  structure(
    list(deltaPhaseRadians = deltaPhaseRadians, ampRatio = ampRatio),
    amplitudes = fitAmp,
    phases     = fitPhase
  )
}

