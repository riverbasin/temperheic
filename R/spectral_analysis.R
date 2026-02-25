# spectral_analysis.R — Power spectrum computation and peak identification
#
# Provides exploratory spectral tools for examining the frequency content
# of temperature time series BEFORE extracting amplitude/phase at specific
# frequencies. Answers the question: "What frequencies carry signal in
# this record?" — a prerequisite for informed use of fit_ols/fit_fft.
#
# Functions:
#   compute_power() — one-sided PSD from a zoo time series (single or multi-column)
#   find_peaks()    — identify dominant spectral peaks by prominence
#
# Design notes:
#   - Returns tidy tibbles (ggplot-ready, pipe-friendly)
#   - Uses base fft() — no new dependencies
#   - Each function is independently useful outside temperheic
#   - Multi-column zoo dispatch: spectra computed per column, returned as
#     long-format tibble with a `series` identifier column


#' Compute power spectral density of a temperature time series
#'
#' Estimates the power spectral density (PSD) of a zoo time series using the
#' periodogram method (base FFT). Returns a tidy tibble for pipe-friendly
#' exploration and ggplot-ready visualization. This is the exploratory
#' complement to \code{\link{fit_ols}} and \code{\link{fit_fft}}: use
#' \code{compute_power()} to survey what frequencies are present in a signal,
#' then use the fitting functions to extract amplitude and phase at specific
#' target frequencies.
#'
#' @details
#' **Method:** The function computes a standard periodogram via \code{stats::fft()}.
#' The input is mean-centered (DC removed) before transforming. An optional
#' Hann window can be applied to reduce spectral leakage from non-integer
#' periodicities, at the cost of ~1.5x broadening of spectral peaks.
#'
#' **Normalization:** Power is normalized as power spectral density (PSD)
#' with units of amplitude\eqn{^2} / Hz, computed as
#' \eqn{2 |X_k|^2 / (N \cdot f_s)} where \eqn{f_s} is the sampling
#' frequency. The factor of 2 accounts for the one-sided spectrum (positive
#' frequencies only). When \code{window = "hann"}, an additional correction
#' factor compensates for the windowing energy loss. PSD values are
#' comparable across records of different lengths and sampling rates.
#'
#' **NA handling:** NAs are filled by linear interpolation
#' (\code{zoo::na.approx}); remaining edge NAs are replaced with the series
#' mean. A warning is issued if more than 10\% of values are NA.
#'
#' **Multi-column zoo objects:** When \code{x} has multiple columns (e.g.,
#' multiple sensor depths), spectra are computed independently for each
#' column and returned as a single long-format tibble with a \code{series}
#' column identifying the source.
#'
#' **Interpreting the output:** The key columns for exploration are
#' \code{period_hr} (or \code{period_day}) and \code{psd}. Plot these on
#' log-log axes to identify dominant periodicities. In temperature data,
#' expect peaks near 24 hr (diel), 12 hr (semi-diurnal), and 8766 hr
#' (annual, if the record is long enough). The relative magnitude of the
#' diel peak across depths is a first-order indicator of signal attenuation
#' and thus subsurface thermal transport.
#'
#' @param x A zoo object. Single-column or multi-column (one per sensor).
#'   Index must be numeric (seconds) — the standard temperheic convention.
#' @param window Character: window function to apply before FFT. One of
#'   \code{"none"} (default, rectangular window) or \code{"hann"} (Hann
#'   window, reduces spectral leakage). The Hann window is recommended when
#'   the record length is not an integer multiple of the target period.
#' @param detrend Logical: if \code{TRUE}, subtract a linear trend before
#'   computing the spectrum. Default \code{FALSE}. Use this when the signal
#'   has a strong non-periodic drift that would leak energy across all
#'   frequencies. For temperature data where you want to see the annual peak,
#'   leave this \code{FALSE}; for examining diel structure in a short window
#'   with residual seasonal curvature, set to \code{TRUE}.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{frequency_hz}{Numeric: frequency in Hz (cycles per second).}
#'     \item{period_s}{Numeric: period in seconds (1 / frequency_hz).}
#'     \item{period_hr}{Numeric: period in hours.}
#'     \item{period_day}{Numeric: period in days.}
#'     \item{psd}{Numeric: power spectral density (amplitude^2 / Hz).}
#'     \item{series}{Character: column name from input zoo object. Present
#'       for all inputs (single-column zoo objects use the column name or
#'       \code{"value"} if unnamed).}
#'   }
#'   Rows are sorted by increasing frequency. Only positive frequencies
#'   are returned (DC and Nyquist components are excluded).
#'
#' @examples
#' \dontrun{
#' # Survey the spectrum of a stream temperature record
#' stream_psd <- stream_zoo %>% compute_power()
#'
#' # Log-log plot of power spectrum
#' library(ggplot2)
#' ggplot(stream_psd, aes(period_hr, psd)) +
#'   geom_line() +
#'   scale_x_log10() + scale_y_log10() +
#'   geom_vline(xintercept = c(12, 24, 8766), lty = 2, alpha = 0.5) +
#'   labs(x = "Period (hours)", y = "PSD")
#'
#' # Compare spectra across sensor depths
#' multi_psd <- multi_sensor_zoo %>% compute_power(window = "hann")
#' ggplot(multi_psd, aes(period_hr, psd, color = series)) +
#'   geom_line() +
#'   scale_x_log10() + scale_y_log10()
#' }
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its
#' Applications (4th ed.). Springer. Chapter 4: Spectral Analysis and Filtering.
#'
#' Luce, C. H., & Tonina, D. (2017). Modeling of streambed surface and
#' hyporheic temperatures. In D. Tonina (Ed.), Open Channel Hydraulics,
#' River Hydraulic Structures and Fluvial Geomorphology. CRC Press.
#' @export
compute_power <- function(x,
                          window = c("none", "hann"),
                          detrend = FALSE) {

  if (!zoo::is.zoo(x)) stop("x must be a zoo object.")
  window <- match.arg(window)

  # Multi-column dispatch: compute per column, bind into long-format tibble
  if (NCOL(x) > 1) {

    col_names <- colnames(x)
    if (is.null(col_names)) col_names <- paste0("V", seq_len(NCOL(x)))

    results <- purrr::imap(
      stats::setNames(as.list(as.data.frame(x)), col_names),
      function(vals, nm) {
        z_single <- zoo::zoo(vals, order.by = zoo::index(x))
        out <- .compute_power_single(z_single, window = window,
                                     detrend = detrend)
        out$series <- nm
        out
      }
    )

    return(dplyr::bind_rows(results))
  }

  # Single-column path
  out <- .compute_power_single(x, window = window, detrend = detrend)
  nm <- colnames(x)
  out$series <- if (!is.null(nm) && length(nm) == 1) nm else "value"
  out
}


#' Compute PSD for a single zoo time series (internal)
#'
#' Workhorse called by \code{compute_power()}. Handles NA interpolation,
#' optional detrending, windowing, FFT, and PSD normalization.
#'
#' @param x A single-column zoo object with numeric (seconds) index.
#' @param window,detrend See \code{compute_power()}.
#' @return A tibble with columns frequency_hz, period_s, period_hr,
#'   period_day, psd. No series column (added by caller).
#' @keywords internal
.compute_power_single <- function(x, window, detrend) {

  idx <- as.numeric(zoo::index(x))
  vals <- as.numeric(zoo::coredata(x))
  n <- length(vals)

  # --- NA handling ---
  n_na <- sum(is.na(vals))
  if (n_na > 0) {
    if (n_na / n > 0.10) {
      warning(sprintf(
        "compute_power: %.1f%% of values are NA. Spectral estimates may be unreliable.",
        100 * n_na / n
      ))
    }
    # Linear interpolation for interior NAs
    x_clean <- zoo::na.approx(x, na.rm = FALSE)
    vals <- as.numeric(zoo::coredata(x_clean))
    # Fill remaining edge NAs with mean
    vals[is.na(vals)] <- mean(vals, na.rm = TRUE)
  }

  # Sampling interval and frequency (Hz)
  dt <- stats::median(diff(idx))
  if (dt <= 0) {
    stop("Non-positive median sampling interval. Check for duplicate timestamps.")
  }
  fs <- 1 / dt

  # --- Preprocessing ---
  # Remove mean (DC component)
  vals <- vals - mean(vals)

  # Optional linear detrend: removes secular drift that would concentrate
  # power at the lowest frequency bins
  if (detrend) {
    t_num <- idx - idx[1]
    fit <- stats::lm.fit(x = cbind(1, t_num), y = vals)
    vals <- fit$residuals
  }

  # --- Window function ---
  # Reduces spectral leakage from non-integer periodicities.
  # Energy correction applied below so PSD amplitude remains unbiased.
  window_correction <- 1.0
  if (window == "hann") {
    w <- 0.5 * (1 - cos(2 * pi * seq(0, n - 1) / (n - 1)))
    vals <- vals * w
    window_correction <- mean(w^2)
  }

  # --- FFT and one-sided PSD ---
  fft_result <- stats::fft(vals)

  # Positive frequencies only, excluding DC (index 1) and Nyquist fold
  nyq <- floor(n / 2)
  pos_idx <- 2:nyq

  # Frequency vector for bins 2:nyq
  freqs <- ((1:(nyq - 1)) * fs) / n

  # One-sided PSD: 2 * |X[k]|^2 / (N * fs)
  # Factor of 2 folds negative-frequency energy into the positive side.
  # Dividing by (N * fs) gives units of amplitude^2 / Hz (standard periodogram).
  psd_raw <- (2 * Mod(fft_result[pos_idx])^2) / (n * fs)

  # Correct for window energy loss
  psd <- psd_raw / window_correction

  # --- Build output tibble ---
  tibble::tibble(
    frequency_hz = freqs,
    period_s     = 1 / freqs,
    period_hr    = 1 / freqs / 3600,
    period_day   = 1 / freqs / 86400,
    psd          = psd
  )
}


#' Identify dominant spectral peaks by prominence
#'
#' Finds local maxima in a power spectrum that exceed a prominence
#' threshold (multiple of the median spectral power). Returns peaks
#' ranked by power, with amplitude computed from the PSD-to-amplitude
#' relationship for a pure cosine.
#'
#' @details
#' **Prominence threshold:** A peak must exceed \code{prominence} times
#' the median PSD to be retained. Higher values yield fewer, more
#' confident peaks. For typical stream temperature spectra,
#' \code{prominence = 100} (default) captures the diel and seasonal
#' peaks while suppressing noise. Increase to 1000+ for very clean
#' signals; decrease to 10-50 for weak or noisy records.
#'
#' **Minimum separation:** When \code{min_separation_hr} is specified,
#' peaks are selected greedily from strongest to weakest, rejecting any
#' candidate within \code{min_separation_hr} hours of an already-selected
#' peak. This prevents multiple adjacent bins around a single true peak
#' from appearing as separate detections.
#'
#' **Amplitude interpretation:** The amplitude column converts PSD back
#' to signal amplitude. For a pure cosine of amplitude A, the PSD peak
#' equals A^2 / (2 * df) where df is the frequency resolution. The
#' amplitude column computes \code{sqrt(psd * df * 2)} accordingly.
#' For broad peaks (energy spread across bins due to non-stationarity
#' or windowing), the reported amplitude is a lower bound.
#'
#' @param power_tbl A tibble from \code{\link{compute_power}()}.
#' @param prominence Numeric: minimum ratio of peak PSD to median PSD.
#'   Default 100.
#' @param min_separation_hr Numeric or NULL: minimum period separation
#'   (hours) between retained peaks. Default NULL (no separation filter).
#' @param min_period_hr Numeric or NULL: discard peaks with period below
#'   this value (hours). Useful for ignoring high-frequency noise.
#' @param max_period_hr Numeric or NULL: discard peaks with period above
#'   this value (hours).
#' @param n_peaks Integer: maximum number of peaks to return. Default 20.
#'
#' @return A tibble with the same columns as \code{compute_power()} output
#'   plus:
#'   \describe{
#'     \item{amplitude}{Estimated amplitude (signal units, e.g. degrees C)
#'       computed from PSD and frequency resolution.}
#'   }
#'   Rows are ordered by decreasing PSD (strongest peak first).
#'
#' @examples
#' \dontrun{
#' # Find dominant frequencies in a stream temperature record
#' spectrum <- stream_zoo %>% compute_power()
#' peaks <- spectrum %>% find_peaks(min_separation_hr = 2)
#' print(peaks)
#'
#' # Focus on the diel band (6-48 hours)
#' diel_peaks <- spectrum %>%
#'   find_peaks(min_period_hr = 6, max_period_hr = 48, prominence = 50)
#' }
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and
#' Its Applications (4th ed.), Chapter 4. Springer.
#' @export
find_peaks <- function(power_tbl,
                       prominence = 100,
                       min_separation_hr = NULL,
                       min_period_hr = NULL,
                       max_period_hr = NULL,
                       n_peaks = 20L) {

  # Validate input: must have the columns produced by compute_power()
  required_cols <- c("frequency_hz", "period_hr", "psd")
  if (!all(required_cols %in% names(power_tbl))) {
    stop("power_tbl must have columns: ",
         paste(required_cols, collapse = ", "),
         ". Use compute_power() to generate it.")
  }

  tbl <- power_tbl
  if (!is.null(min_period_hr)) tbl <- tbl[tbl$period_hr >= min_period_hr, ]
  if (!is.null(max_period_hr)) tbl <- tbl[tbl$period_hr <= max_period_hr, ]

  pwr <- tbl$psd
  n <- length(pwr)

  if (n < 3) {
    out <- tbl[integer(0), ]
    out$amplitude <- numeric(0)
    return(out)
  }

  # --- Local maxima: psd[i] > both neighbors ---
  is_peak <- rep(FALSE, n)
  is_peak[2:(n - 1)] <- (pwr[2:(n - 1)] > pwr[1:(n - 2)]) &
                         (pwr[2:(n - 1)] > pwr[3:n])

  # --- Prominence filter ---
  median_psd <- stats::median(pwr, na.rm = TRUE)
  is_prominent <- pwr > (prominence * median_psd)

  candidates <- which(is_peak & is_prominent)
  if (length(candidates) == 0) {
    out <- tbl[integer(0), ]
    out$amplitude <- numeric(0)
    return(out)
  }

  # Sort by PSD, strongest first (greedy selection)
  candidates <- candidates[order(-pwr[candidates])]

  # --- Minimum separation filter ---
  if (!is.null(min_separation_hr)) {
    periods <- tbl$period_hr
    kept <- integer(0)
    for (ci in candidates) {
      if (length(kept) == 0 ||
          all(abs(periods[ci] - periods[kept]) > min_separation_hr)) {
        kept <- c(kept, ci)
      }
      if (length(kept) >= n_peaks) break
    }
    candidates <- kept
  } else {
    candidates <- utils::head(candidates, n_peaks)
  }

  out <- tbl[candidates, ]

  # Amplitude from PSD: for a pure cosine of amplitude A,

  # PSD_peak = A^2 / (2 * df) where df = 1/(N*dt) is frequency resolution.
  # Since we already divided by df in the PSD calculation,
  # amplitude = sqrt(psd * df * 2) is the inverse.
  # However, for the user the simpler approximation using the frequency
  # resolution from the data is more transparent:
  df <- tbl$frequency_hz[2] - tbl$frequency_hz[1]
  out$amplitude <- sqrt(out$psd * df * 2)

  out[order(-out$psd), ]
}

