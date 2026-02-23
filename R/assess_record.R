# assess_record.R — Evaluate fitness of a temperature record for periodic fitting
#
# Diagnoses whether a zoo time series has sufficient length and continuity
# to support reliable amplitude/phase extraction at a target period. Based
# on Johnson et al. (2021) "Heed the data gap" empirical thresholds for
# annual OLS sine-wave fits, generalized to arbitrary target periods via
# a frequency-domain framing.
#
# For diel analysis, gap tolerance and segment sufficiency interact with
# the windowed fitting workflow (Vogt et al. 2012, Gariglio et al. 2013).
# The min_segment_duration parameter allows users to enforce a minimum
# usable segment length independent of the cycle-count criterion.


#' Assess whether a temperature record supports fitting at a target period
#'
#' Examines a zoo time series for gaps and total length relative to a target
#' period (diel, annual, or arbitrary). Returns a diagnostic tibble describing
#' contiguous segments, gap locations and durations, and a per-segment verdict
#' on whether the record is sufficient for OLS or FFT fitting.
#'
#' @details
#' The assessment applies two period-relative criteria, plus an optional
#' absolute duration floor:
#'
#' **Maximum gap tolerance** (\code{max_gap_frac}): The longest permissible gap as
#' a fraction of the target period. Gaps shorter than this are tolerable
#' (interpolation is defensible); gaps longer than this split the record into
#' separate contiguous segments. Default 0.17 (\eqn{\approx} 62 days for annual),
#' based on Johnson et al. (2021) finding that up to 7--9 consecutive weeks of
#' missing data still yields accurate annual sine-wave parameter estimates.
#'
#' **Minimum record length** (\code{min_cycles}): The minimum number of complete
#' target-period cycles a contiguous segment must span. Segments shorter than
#' \code{min_cycles * period} are flagged as insufficient.
#'
#' **Minimum segment duration** (\code{min_segment_duration}): An optional absolute
#' floor on segment duration in seconds, independent of cycle count. When set,
#' a segment must satisfy \emph{both} \code{min_cycles} and \code{min_segment_duration}
#' to be classified as "sufficient". This is particularly useful for diel
#' analysis, where \code{min_cycles = 5} yields a 5-day minimum, but the user may
#' want to enforce a longer minimum to ensure enough data for windowed fitting
#' (e.g., \code{min_segment_duration = 7 * 86400} for 7-day windows).
#'
#' Default thresholds and their basis:
#' \itemize{
#'   \item \code{max_gap_frac = 0.17}: 62 days / 365 days. Johnson et al. (2021)
#'     showed OLS sine fits remain accurate with up to ~9 weeks missing
#'     (without imputation). With imputation, tolerance extends to ~13 weeks
#'     (0.25 of period).
#'   \item \code{min_cycles = 1.5} (OLS minimum): Below 1.5 cycles, the seasonal
#'     curvature is confounded with trend and amplitude is poorly constrained.
#'   \item \code{min_cycles = 2.0} (FFT minimum): The Rayleigh criterion requires
#'     \eqn{\geq 2} cycles to resolve a frequency as a distinct spectral peak.
#'   \item \code{min_cycles = 3.0} (recommended): Empirical rule of thumb for
#'     reliable amplitude estimates (\eqn{\pm \sim 10\%}).
#' }
#'
#' **Gap detection:** A gap is any interval between consecutive observations
#' that exceeds \code{gap_threshold} times the median sampling interval. The default
#' multiplier of 3 catches instrument outages while tolerating occasional
#' missing hourly records.
#'
#' **Diel analysis considerations:** For windowed diel analysis (e.g., using
#' \code{window_temperature()} with 5--7 day windows per Vogt et al. 2012 and
#' Gariglio et al. 2013), the practical minimum segment length is governed by
#' the window width, not the cycle count alone. A segment of 3 diel cycles
#' (3 days) satisfies \code{min_cycles = 2.0} but is too short for a 5-day analysis
#' window. Use \code{min_segment_duration} to enforce the window-width floor:
#'
#' \preformatted{
#'   # Diel assessment requiring 7-day minimum segments for windowed fitting
#'   assess_record(z, period = 86400, min_cycles = 5.0,
#'                 min_segment_duration = 7 * 86400)
#' }
#'
#' The default \code{max_gap_frac = 0.17} translates to \eqn{\approx 4.1} hours
#' for diel period. This is a reasonable tolerance: a gap shorter than
#' ~4 hours leaves most of the diel sinusoid constrained for OLS fitting
#' within each window. Longer gaps (e.g., 10+ hours) lose substantial
#' curvature information and should split the record. Note that
#' \code{window_temperature()} has its own \code{min_coverage} parameter that
#' provides within-window gap checking; \code{assess_record()} provides the
#' complementary record-level diagnostic.
#'
#' @param x A zoo object (single-column). Index should be numeric (seconds)
#'   or POSIXct.
#' @param period Numeric: target period in seconds. Common values:
#'   86400 (diel), 31557600 (annual = 365.25 * 86400). No default; user
#'   must specify.
#' @param max_gap_frac Numeric: maximum tolerable gap as a fraction of
#'   \code{period}. Default 0.17 (~9 weeks for annual period, per Johnson et al.
#'   2021; ~4.1 hours for diel). Set to 0.25 if imputation will be applied.
#' @param min_cycles Numeric: minimum contiguous record length expressed as
#'   multiples of \code{period}. Default 2.0 (Nyquist/Rayleigh minimum for FFT).
#'   Use 1.5 for OLS-only, 3.0 for confident estimates. For diel windowed
#'   analysis, consider 5.0+ to ensure enough cycles for a full analysis
#'   window.
#' @param gap_threshold Numeric: multiplier on median sampling interval to
#'   detect gaps. Default 3. A gap is flagged when consecutive observations
#'   are separated by more than \code{gap_threshold * median(diff(index))}.
#' @param min_segment_duration Numeric or NULL: optional absolute minimum
#'   segment duration in seconds. When non-NULL, a segment must meet both
#'   \code{min_cycles * period} and \code{min_segment_duration} to be classified as
#'   "sufficient". Useful for diel analysis where windowed fitting requires
#'   a minimum number of days regardless of cycle count. Default NULL
#'   (no absolute floor; classification depends only on \code{min_cycles}).
#'
#' @return A list with class \code{"record_assessment"} containing:
#'   \describe{
#'     \item{segments}{A tibble with one row per contiguous segment:
#'       \code{seg_id}, \code{start}, \code{end}, \code{duration_s}, \code{n_obs}, \code{cycles}
#'       (duration / period), \code{verdict} ("sufficient", "marginal",
#'       "insufficient").}
#'     \item{gaps}{A tibble with one row per detected gap:
#'       \code{gap_id}, \code{start}, \code{end}, \code{duration_s}, \code{frac_of_period},
#'       \code{tolerable} (logical: duration < max_gap_frac * period).}
#'     \item{overall_verdict}{Character: "sufficient" if any segment meets
#'       \code{min_cycles} (and \code{min_segment_duration} if set), "marginal" if
#'       any segment >= 1.0 cycles but none meets full sufficiency criteria,
#'       "insufficient" otherwise.}
#'     \item{params}{List of the parameter values used: \code{period},
#'       \code{max_gap_frac}, \code{min_cycles}, \code{gap_threshold}, \code{max_gap_s},
#'       \code{min_record_s}, \code{min_segment_duration}.}
#'   }
#'
#' @references
#' Johnson, Z. C., Johnson, B. G., Briggs, M. A., Snyder, C. D., Hitt,
#' N. P., & Devine, W. D. (2021). Heed the data gap: Guidelines for using
#' incomplete datasets in annual stream temperature analyses. Ecological
#' Indicators, 122, 107229. doi:10.1016/j.ecolind.2020.107229
#'
#' Hare, D. K., Helton, A. M., Johnson, Z. C., Lane, J. W., & Briggs,
#' M. A. (2023). PASTA: Practical Automated Stream Temperature Analysis.
#' Water Resources Research, 59, e2022WR033912. doi:10.1029/2022WR033912
#'
#' Vogt, T., Schirmer, M., & Cirpka, O. A. (2012). Investigating riparian
#' groundwater flow close to a losing river using diurnal temperature
#' oscillations at high vertical resolution. Hydrology and Earth System
#' Sciences, 16, 473--487. doi:10.5194/hess-16-473-2012
#'
#' Gariglio, F. P., Tonina, D., & Luce, C. H. (2013). Spatiotemporal
#' variability of hyporheic exchange through a pool-riffle-pool sequence.
#' Water Resources Research, 49, 7185--7204. doi:10.1002/wrcr.20419
#'
#' @examples
#' \dontrun{
#' # Assess a spring temperature record for annual fitting
#' assessment <- assess_record(spring6_full_zoo, period = 365.25 * 86400)
#' print(assessment)
#'
#' # Stricter threshold for FFT
#' assess_record(well10_full_zoo, period = 365.25 * 86400, min_cycles = 3.0)
#'
#' # Diel assessment with default gap tolerance
#' assess_record(stream_zoo, period = 86400, min_cycles = 5.0)
#'
#' # Diel assessment enforcing 7-day minimum segments for windowed analysis
#' assess_record(stream_zoo, period = 86400, min_cycles = 5.0,
#'               min_segment_duration = 7 * 86400)
#' }
#'
#' @export
assess_record <- function(x,
                          period,
                          max_gap_frac = 0.17,
                          min_cycles = 2.0,
                          gap_threshold = 3,
                          min_segment_duration = NULL) {

  if (!zoo::is.zoo(x)) stop("x must be a zoo object.")
  if (missing(period)) stop("period (seconds) must be specified.")
  if (period <= 0) stop("period must be positive.")
  if (!is.null(min_segment_duration) && min_segment_duration <= 0) {
    stop("min_segment_duration must be positive (seconds) or NULL.")
  }
  if (NCOL(x) > 1) {
    warning("Multi-column zoo: assessing based on index only (gaps are shared).",
            " Columns are not evaluated independently.")
  }

  # --- Derived thresholds ---
  max_gap_s    <- max_gap_frac * period
  min_record_s <- min_cycles * period

  # Effective minimum segment duration is the more restrictive of the two
  # criteria: cycle-based and absolute duration floor (if set).
  # A segment must exceed BOTH to qualify as "sufficient".
  effective_min_s <- if (!is.null(min_segment_duration)) {
    max(min_record_s, min_segment_duration)
  } else {
    min_record_s
  }

  # --- Index handling ---
  idx <- as.numeric(zoo::index(x))
  n   <- length(idx)

  if (n < 2) {
    stop("Record has fewer than 2 observations.")
  }

  # --- Gap detection ---
  # A gap is any interval > gap_threshold * median sampling interval.
  # This catches instrument outages while tolerating occasional missing
  # records (e.g., 1-2 missing hourly obs in an otherwise continuous record).
  dt       <- diff(idx)
  median_dt <- stats::median(dt)
  gap_locs <- which(dt > gap_threshold * median_dt)

  # --- Build segment and gap tables ---
  # Segments are bounded by gaps (or record start/end)
  seg_starts <- c(1, gap_locs + 1)
  seg_ends   <- c(gap_locs, n)
  n_segs     <- length(seg_starts)

  segments <- tibble::tibble(
    seg_id     = seq_len(n_segs),
    start      = idx[seg_starts],
    end        = idx[seg_ends],
    duration_s = idx[seg_ends] - idx[seg_starts],
    n_obs      = seg_ends - seg_starts + 1L,
    cycles     = (idx[seg_ends] - idx[seg_starts]) / period
  )

  # Verdict per segment: must satisfy both cycle count AND absolute
  # duration floor (if set) to be "sufficient"
  segments$verdict <- dplyr::case_when(
    segments$cycles >= min_cycles &
      segments$duration_s >= effective_min_s ~ "sufficient",
    segments$cycles >= 1.0                   ~ "marginal",
    TRUE                                     ~ "insufficient"
  )

  # Gap table
  if (length(gap_locs) > 0) {
    gaps <- tibble::tibble(
      gap_id          = seq_along(gap_locs),
      start           = idx[gap_locs],
      end             = idx[gap_locs + 1],
      duration_s      = idx[gap_locs + 1] - idx[gap_locs],
      frac_of_period  = (idx[gap_locs + 1] - idx[gap_locs]) / period,
      tolerable       = (idx[gap_locs + 1] - idx[gap_locs]) < max_gap_s
    )
  } else {
    gaps <- tibble::tibble(
      gap_id         = integer(0),
      start          = numeric(0),
      end            = numeric(0),
      duration_s     = numeric(0),
      frac_of_period = numeric(0),
      tolerable      = logical(0)
    )
  }

  # --- Overall verdict ---
  if (any(segments$verdict == "sufficient")) {
    overall <- "sufficient"
  } else if (any(segments$verdict == "marginal")) {
    overall <- "marginal"
  } else {
    overall <- "insufficient"
  }

  structure(
    list(
      segments        = segments,
      gaps            = gaps,
      overall_verdict = overall,
      params          = list(
        period               = period,
        max_gap_frac         = max_gap_frac,
        min_cycles           = min_cycles,
        gap_threshold        = gap_threshold,
        max_gap_s            = max_gap_s,
        min_record_s         = min_record_s,
        min_segment_duration = min_segment_duration
      )
    ),
    class = "record_assessment"
  )
}


#' Print method for record_assessment objects
#'
#' @param x A record_assessment object from \code{assess_record()}.
#' @param ... Additional arguments (ignored).
#' @export
print.record_assessment <- function(x, ...) {

  p <- x$params

  # Human-readable period name
  period_days <- p$period / 86400
  period_label <- if (abs(period_days - 1) < 0.01) {
    "diel (1 day)"
  } else if (abs(period_days - 365.25) < 1) {
    "annual (365.25 days)"
  } else {
    paste0(round(period_days, 1), " days")
  }

  cat("Record assessment for", period_label, "period
")
  cat("  Max gap tolerance:", round(p$max_gap_s / 86400, 1), "days",
      paste0("(", round(p$max_gap_frac * 100, 0), "% of period)
"))
  cat("  Min record length:", round(p$min_record_s / 86400, 0), "days",
      paste0("(", p$min_cycles, " cycles)
"))
  if (!is.null(p$min_segment_duration)) {
    cat("  Min segment duration:", round(p$min_segment_duration / 86400, 1),
        "days (absolute floor)
")
  }
  cat("
")

  # Segments
  cat("Segments (", nrow(x$segments), "):
", sep = "")
  for (i in seq_len(nrow(x$segments))) {
    s <- x$segments[i, ]
    t_start <- format(as.POSIXct(s$start, origin = "1970-01-01", tz = "UTC"),
                      "%Y-%m-%d")
    t_end   <- format(as.POSIXct(s$end, origin = "1970-01-01", tz = "UTC"),
                      "%Y-%m-%d")
    cat(sprintf("  [%d] %s to %s  %.0f days  %.2f cycles  [%s]
",
                s$seg_id, t_start, t_end,
                s$duration_s / 86400, s$cycles,
                toupper(s$verdict)))
  }

  # Gaps
  if (nrow(x$gaps) > 0) {
    cat("
Gaps (", nrow(x$gaps), "):
", sep = "")
    for (i in seq_len(nrow(x$gaps))) {
      g <- x$gaps[i, ]
      t_start <- format(as.POSIXct(g$start, origin = "1970-01-01", tz = "UTC"),
                        "%Y-%m-%d")
      t_end   <- format(as.POSIXct(g$end, origin = "1970-01-01", tz = "UTC"),
                        "%Y-%m-%d")
      flag <- if (g$tolerable) "ok" else "EXCEEDS TOLERANCE"
      cat(sprintf("  [%d] %s to %s  %.1f days  (%.0f%% of period)  [%s]
",
                  g$gap_id, t_start, t_end,
                  g$duration_s / 86400,
                  g$frac_of_period * 100,
                  flag))
    }
  } else {
    cat("
No gaps detected.
")
  }

  cat("
Overall verdict:", toupper(x$overall_verdict), "
")
  invisible(x)
}

