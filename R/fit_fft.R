
#' Fit amplitude and phase using FFT spectral extraction
#'
#' Extracts amplitude and phase at the nearest Fourier frequency bin to the
#' target period. A single FFT gives amplitude and phase at every resolvable
#' frequency simultaneously; this function extracts the bin closest to the
#' specified period.
#'
#' @details
#' **Phase convention:** To match fit_ols, phase is computed as -Arg(X[k])
#' where X[k] is the complex DFT coefficient at the target bin. This gives
#' phase in the same convention as atan2(alpha_s, alpha_c) from the OLS
#' linearized cosine model, ensuring pairwise deltaPhaseRadians are
#' consistent across methods.
#'
#' **Bin matching:** The target period (e.g., 86400s for diel) may not land
#' exactly on a Fourier frequency bin. The function selects the nearest bin.
#' Accuracy improves with longer records (finer frequency resolution).
#' A warning is issued if the nearest bin period differs from the target
#' by more than 1%.
#'
#' **When to prefer FFT over OLS:**
#' - When you want to survey multiple frequencies from the same data
#' - When the signal contains harmonics (non-sinusoidal shape) and you
#'   want amplitude/phase at each harmonic independently
#' - As a cross-check against OLS results
#'
#' **When to prefer OLS:**
#' - When the target frequency does not land near a Fourier bin
#' - For short records where frequency resolution is coarse
#' - When you need to include polynomial trend terms in the regression
#'
#' @param empiricalData A zoo object with one column per sensor depth.
#' @param boundaryMean Numeric -- not used by FFT but retained for interface
#'   compatibility.
#' @param periodInSeconds Numeric -- period of the target cycle (s).
#' @param optimizeRange Numeric length-2 -- phase window for unwrapping.
#' @param nmin Integer -- minimum non-NA observations to attempt fitting.
#' @param empiricalDataPeriods Numeric -- number of complete cycles.
#'
#' @return Same structure as fit_ols/fitCosine: list(deltaPhaseRadians,
#'   ampRatio) with attributes "amplitudes" and "phases".
#'
#' @references
#' Luce, C. H., Tonina, D., Gariglio, F., & Applebee, R. (2013). WRR 49.
#' Vogt, T., Schirmer, M., & Cirpka, O. A. (2012). HESS 16.
#' @export
fit_fft <- function(empiricalData, boundaryMean, periodInSeconds,
                    optimizeRange, nmin, empiricalDataPeriods) {

  # Time vector in seconds from start
  seconds <- as.numeric(zoo::index(empiricalData))
  obsTime <- seconds - seconds[1]
  n <- length(obsTime)
  dt <- stats::median(diff(obsTime))

  # Fourier frequency vector (Hz) and target bin lookup
  freqs <- (0:(n - 1)) / (n * dt)
  target_freq <- 1 / periodInSeconds
  nyquist <- floor(n / 2)
  target_bin <- which.min(abs(freqs[1:nyquist] - target_freq))

  # Warn if bin mismatch exceeds 1%
  actual_period <- 1 / freqs[target_bin]
  bin_mismatch <- abs(actual_period - periodInSeconds) / periodInSeconds
  if (bin_mismatch > 0.01) {
    warning(sprintf(
      paste0("fit_fft: nearest FFT bin period (%.1f s) differs from target ",
             "(%.1f s) by %.1f%%. Consider a longer record for better ",
             "frequency resolution."),
      actual_period, periodInSeconds, bin_mismatch * 100
    ))
  }

  # --- Extract amplitude and phase at target bin for each sensor ---
  fits <- purrr::map(as.data.frame(empiricalData), function(eD) {

    # Require minimum non-NA observations
    non_na <- stats::na.omit(eD)
    if (length(non_na) < nmin) {
      return(c(fitAmp = NA_real_, fitPhase = NA_real_))
    }

    # Interpolate NAs for FFT (requires complete vector)
    eD_clean <- zoo::na.approx(eD, na.rm = FALSE)
    # Fill any remaining edge NAs with mean
    eD_clean[is.na(eD_clean)] <- mean(eD_clean, na.rm = TRUE)

    # FFT
    fft_result <- stats::fft(eD_clean)
    coeff <- fft_result[target_bin]

    # Amplitude: 2 * |X[k]| / N (factor of 2 because energy split pos/neg freq)
    amp <- 2 * Mod(coeff) / n

    # Phase: -Arg(X[k]) to match OLS convention (atan2(alpha_s, alpha_c))
    pha <- (-Arg(coeff)) %% (2 * pi)

    c(fitAmp = amp, fitPhase = pha)
  })

  # Extract per-sensor amplitude and phase vectors
  fitAmp   <- purrr::map_dbl(fits, "fitAmp")
  fitPhase <- purrr::map_dbl(fits, "fitPhase")

  # --- Phase unwrapping (same logic as fit_ols) ---
  initialPhase <- fitPhase[1]
  relativeRange <- 2 * pi * optimizeRange + initialPhase
  fitPhase[fitPhase < relativeRange[1]] <- fitPhase[fitPhase < relativeRange[1]] + 2 * pi
  fitPhase <- fitPhase + (empiricalDataPeriods - 1) * 2 * pi

  # --- Pairwise amplitude ratios and phase differences ---
  combos <- expand.grid(from = seq_along(fitPhase), to = seq_along(fitAmp))
  deltaPhaseRadians <- fitPhase[combos$to] - fitPhase[combos$from]
  ampRatio          <- fitAmp[combos$to]   / fitAmp[combos$from]

  structure(
    list(deltaPhaseRadians = deltaPhaseRadians, ampRatio = ampRatio),
    amplitudes = fitAmp,
    phases     = fitPhase
  )
}

