# signal_processing.R — Preprocessing functions for temperature time series
#
# Provides atomic, composable tools for preparing field temperature data
# before amplitude/phase extraction. Each function accepts and returns zoo
# objects, making them independently useful outside the temperheic pipeline.
#
# Functions:
#   detrend_temperature() — remove slow trend to isolate oscillatory component


#' Remove trend from a temperature time series
#'
#' Estimates and subtracts a slowly-varying trend component, isolating the
#' oscillatory signal (diel, seasonal, or other periodic component) for
#' subsequent amplitude/phase extraction. The trend is typically the annual
#' cycle when targeting diel signals, or a multi-year drift when targeting
#' seasonal signals.
#'
#' @details
#' Three detrending methods are available:
#'
#' **loess** (default): Local polynomial regression via `stats::loess()`.
#'   The `span` parameter controls the smoothing bandwidth as a fraction of
#'   the data. Smaller span = more flexible trend. A span of ~0.05-0.1
#'   captures seasonal variation in a year-long record while preserving
#'   diel oscillations. This is the most flexible method and generally
#'   the best default for field data with irregular trends.
#'
#' **ma** (moving average): Centered running mean via `stats::filter()`.
#'   The `window` parameter specifies the width in number of observations.
#'   For hourly data targeting the diel cycle, a 24-point window (1 day)
#'   removes the diel signal and leaves the trend; a 168-point window
#'   (1 week) smooths both diel and synoptic variability. Edge values are
#'   NA; use `na.rm = TRUE` in downstream functions.
#'
#' **polynomial**: Polynomial regression via `stats::lm()`. The `degree`
#'   parameter controls the polynomial order (1 = linear, 2 = quadratic,
#'   etc.). Best for records where the trend is approximately monotonic
#'   or gently curved. Avoids the edge effects of loess/ma but cannot
#'   capture complex seasonal patterns.
#'
#' For multi-column zoo objects (multiple sensor depths), detrending is
#' applied independently to each column via `purrr::map()`.
#'
#' @param x A zoo object. Single-column or multi-column (one per sensor).
#' @param method Character: detrending method. One of "loess" (default),
#'   "ma" (moving average), or "polynomial".
#' @param span Numeric: loess smoothing span as fraction of data (0, 1].
#'   Default 0.05. Only used when `method = "loess"`.
#' @param window Integer: moving average window width in observations.
#'   Default 24 (one day at hourly sampling). Only used when `method = "ma"`.
#' @param degree Integer: polynomial degree. Default 2. Only used when
#'   `method = "polynomial"`.
#' @param return_trend Logical: if TRUE, return a list with both the
#'   detrended series and the estimated trend. Default FALSE (returns
#'   only the detrended zoo object for pipe-friendly use).
#'
#' @return If `return_trend = FALSE` (default): a zoo object of the same
#'   structure as `x`, containing the detrended (residual) series.
#'   If `return_trend = TRUE`: a list with components `detrended` (zoo)
#'   and `trend` (zoo), both matching the structure of `x`.
#'
#' @examples
#' \dontrun{
#' # Remove seasonal trend from a year of hourly stream temperature
#' stream_detrended <- stream_temp %>%
#'   detrend_temperature(method = "loess", span = 0.05)
#'
#' # Inspect both trend and residual
#' result <- stream_temp %>%
#'   detrend_temperature(method = "ma", window = 168, return_trend = TRUE)
#' plot(result$trend)
#' plot(result$detrended)
#' }
#'
#' @export
detrend_temperature <- function(x,
                                method = c("loess", "ma", "polynomial"),
                                span = 0.05,
                                window = 24L,
                                degree = 2L,
                                return_trend = FALSE) {

  method <- match.arg(method)

  if (!zoo::is.zoo(x)) stop("x must be a zoo object.")

  # Dispatch: apply detrending column-wise for multi-column zoo objects,
  # or directly for single-column. purrr::map preserves names.
  if (NCOL(x) > 1) {
    # Multi-column: detrend each sensor independently
    results <- purrr::map(
      as.list(as.data.frame(x)),
      ~ .detrend_single(
        zoo::zoo(.x, order.by = zoo::index(x)),
        method = method, span = span,
        window = window, degree = degree
      )
    )

    if (return_trend) {
      list(
        detrended = do.call(zoo::merge.zoo, purrr::map(results, "detrended")),
        trend     = do.call(zoo::merge.zoo, purrr::map(results, "trend"))
      )
    } else {
      do.call(zoo::merge.zoo, purrr::map(results, "detrended"))
    }

  } else {
    result <- .detrend_single(x, method = method, span = span,
                              window = window, degree = degree)
    if (return_trend) result else result$detrended
  }
}


#' Detrend a single zoo time series (internal)
#'
#' Workhorse function called by `detrend_temperature()`. Always returns
#' a list with both `detrended` and `trend` components.
#'
#' @param x A single-column zoo object.
#' @param method,span,window,degree See `detrend_temperature()`.
#' @return List with `detrended` (zoo) and `trend` (zoo).
#' @keywords internal
.detrend_single <- function(x, method, span, window, degree) {

  idx <- zoo::index(x)
  vals <- zoo::coredata(x)

  # Numeric time for regression — seconds from start
  t_num <- as.numeric(idx) - as.numeric(idx[1])

  trend_vals <- switch(method,

    loess = {
      # loess needs a data frame; predict back onto original time points
      fit <- stats::loess(vals ~ t_num, span = span,
                          na.action = stats::na.exclude)
      stats::predict(fit, newdata = data.frame(t_num = t_num))
    },

    ma = {
      # Centered moving average: sides = 2 for symmetric (no phase shift)
      as.numeric(stats::filter(vals, rep(1 / window, window), sides = 2))
    },

    polynomial = {
      # Polynomial fit: orthogonal polynomials for numerical stability
      fit <- stats::lm(vals ~ stats::poly(t_num, degree = degree, raw = FALSE),
                       na.action = stats::na.exclude)
      stats::predict(fit, newdata = data.frame(t_num = t_num))
    }
  )

  trend <- zoo::zoo(trend_vals, order.by = idx)
  detrended <- zoo::zoo(vals - trend_vals, order.by = idx)

  list(detrended = detrended, trend = trend)
}



#' Slice a temperature time series into overlapping analysis windows
#'
#' Divides a zoo object into temporal windows for time-varying analysis.
#' Each window becomes an independent input for amplitude/phase extraction
#' and flux estimation, producing a time series of parameter estimates
#' rather than a single whole-record average.
#'
#' @details
#' Windows are defined by width and step, both in seconds (matching the
#' zoo index units and the physical constants in the ADE solution).
#' Common values:
#'
#' \itemize{
#'   \item Diel analysis: `width = 7 * 86400` (7 days), `step = 86400` (1 day)
#'   \item Seasonal analysis: `width = 90 * 86400` (90 days), `step = 30 * 86400` (30 days)
#' }
#'
#' The \pkg{lubridate} package provides a readable bridge:
#' `as.numeric(lubridate::as.duration(lubridate::days(7)), "seconds")` returns
#' 604800 (= 7 * 86400). Internalizing these conversions builds intuition
#' for the connection between the time axis, angular frequency (omega =
#' 2*pi/period), and the physical process.
#'
#' Key periods in seconds:
#' \itemize{
#'   \item 1 hour = 3600
#'   \item 1 day (diel period) = 86400
#'   \item 1 week = 604800
#'   \item 1 year (annual period) = 31557600 (365.25 * 86400)
#' }
#'
#' Windows that extend beyond the data are silently dropped. A warning
#' is issued if fewer than 2 complete windows are produced.
#'
#' For multi-column zoo objects, the entire object is windowed together
#' (all sensors share the same time axis), preserving the column structure
#' needed by `thObservedSeries()`.
#'
#' @param x A zoo object. Single-column or multi-column.
#' @param width Numeric: window width in seconds. For diel analysis,
#'   7 * 86400 (7 days) is a good starting point.
#' @param step Numeric: step size in seconds between successive window
#'   starts. Controls overlap. `step = width` gives non-overlapping
#'   windows; `step < width` gives overlap; `step = 86400` slides by
#'   one day. Default is `width` (non-overlapping).
#' @param min_coverage Numeric in (0, 1]: minimum fraction of non-NA
#'   values required within a window to retain it. Default 0.9.
#'   Windows with excessive missing data are dropped with a message.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{window_start}{Numeric: start time of window (seconds, matching zoo index).}
#'     \item{window_end}{Numeric: end time of window (seconds).}
#'     \item{window_center}{Numeric: midpoint of window (seconds).}
#'     \item{data}{List-column of zoo objects, one per window. Each zoo
#'       retains the column structure of the input.}
#'   }
#'
#' @examples
#' \dontrun{
#' # 7-day sliding windows with 1-day step for diel analysis
#' windows <- stream_temp %>%
#'   window_temperature(width = 7 * 86400, step = 86400)
#'
#' # Pipe into analysis with purrr
#' results <- windows %>%
#'   dplyr::mutate(
#'     detrended = purrr::map(data, detrend_temperature),
#'     flux = purrr::map(detrended, ~estimate_flux(.x, ...))
#'   )
#' }
#'
#' @export
window_temperature <- function(x,
                               width,
                               step = width,
                               min_coverage = 0.9) {

  if (!zoo::is.zoo(x)) stop("x must be a zoo object.")
  if (width <= 0) stop("width must be positive.")
  if (step <= 0) stop("step must be positive.")
  if (min_coverage <= 0 || min_coverage > 1) {
    stop("min_coverage must be in (0, 1].")
  }

  idx <- as.numeric(zoo::index(x))
  t_start <- idx[1]
  t_end   <- idx[length(idx)]

  # Generate window start times: first window starts at data start,
  # last window must fit entirely within data range
  if ((t_end - t_start) < width) {
    stop("No complete windows fit within the data. width (", width,
         " s) exceeds data range (", round(t_end - t_start), " s).")
  }
  starts <- seq(from = t_start, to = t_end - width, by = step)

  if (length(starts) < 2) {
    warning("Fewer than 2 complete windows. Consider a smaller width or step.")
  }


  # Build each window: subset zoo by logical index on time range
  windows <- purrr::map(starts, function(s) {
    x[idx >= s & idx < (s + width)]
  })

  # Check coverage: fraction of non-NA values in each window
  coverage <- purrr::map_dbl(windows, function(w) {
    vals <- zoo::coredata(w)
    sum(!is.na(vals)) / length(vals)
  })

  # Filter by min_coverage
  keep <- coverage >= min_coverage
  n_dropped <- sum(!keep)
  if (n_dropped > 0) {
    message(n_dropped, " window(s) dropped due to insufficient data coverage ",
            "(< ", min_coverage * 100, "%).")
  }

  tibble::tibble(
    window_start  = starts[keep],
    window_end    = starts[keep] + width,
    window_center = starts[keep] + width / 2,
    data          = windows[keep]
  )
}
