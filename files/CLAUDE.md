# temperheic — Temperature Signal Analysis for Hyporheic Exchange

An R package implementing analytical methods (Luce et al. 2013, 2017; van Kampen et al. 2022) for estimating groundwater–surface water exchange from temperature time series. Designed for restoration practitioners, consulting hydrologists, and fisheries biologists — not just academic researchers.

## Project Status

Current state: Core inverse engine validated (156 tests, 0 failures), spectral tools working, 3 demo scripts. See @docs/goals-checklist.md for the full phased plan.

Active work: Phase 0 (commit baseline + code review), then Phase 1 (uncertainty quantification).

## Key References

See @docs/refs/luce-2013-reference.md for the core η-method equations
See @docs/refs/luce-2017-reference.md for boundary condition independence
See @docs/refs/van-kampen-2022-reference.md for the LPMLEn frequency domain method
See @docs/refs/bertagnoli-2024-reference.md for iFLOW framework and test cases

Full PDFs are in `refs/` — read them when exact derivations, figures, or tables are needed.

See @docs/design-philosophy.md for the overarching package design principles
See @docs/tidyverse-fp-guidelines.md for R coding conventions
See @docs/methods-review.md for the complete theoretical framework
See @docs/architecture.md for the broader R toolkit vision (Layer 1–5)
See @docs/project-instructions.md for ReacTran/forward model patterns

## Build & Test Commands

```bash
# Load package in development mode
Rscript -e "devtools::load_all()"

# Run full test suite
Rscript -e "devtools::test()"

# Run R CMD check
Rscript -e "devtools::check()"

# Rebuild roxygen documentation
Rscript -e "devtools::document()"

# Install locally
Rscript -e "devtools::install_local(force = TRUE)"
```

## Package Structure

```
R/
  s3Objects.R         # S3 classes: thUnits, thAquifer, thBoundary, thHydro, thSignal
  fit_ols.R           # OLS amplitude/phase extraction (primary fitting engine)
  fit_fft.R           # FFT-based amplitude/phase extraction (cross-check)
  thSeries.R          # Forward model (thSeries) + inverse solution (thObservedSeries)
  thUtils.R           # Legacy utilities: fitCosine, derived2DArray
  signal_processing.R # detrend_temperature, window_temperature
  spectral_analysis.R # compute_power, find_peaks
  assess_record.R     # Data quality assessment for temperature records
  th_test.R           # Round-trip validation: generate → fit → recover parameters
tests/
  testthat/           # 156 tests, organized by function
inst/
  extdata/            # Meacham Creek datasets, SNOTEL air temp
```

## Critical Conventions

### Units — SI base units everywhere
- `hydCond` in m/s (NOT m/d) — passing m/d is the #1 user error
- `period` in seconds: 86400 (diel), 31557600 (annual)
- `tVals` and zoo indices in seconds
- Convert: `K_m_s = K_m_d / 86400`

### Sign conventions
- Positive z: downward into streambed
- Positive v_t: downwelling
- η < 1: downwelling; η = 1: no flow; η > 1: upwelling

### Internal data flow
- zoo objects internally (matrix + time index, efficient for computation)
- tibble for user-facing outputs
- Convert at boundaries between internal and external interfaces

### Fitting engine priority
- OLS (fit_ols): primary — machine-precision recovery, no starting values, closed-form
- FFT (fit_fft): cross-check and spectral exploration
- NLS (fitCosine in thUtils): legacy, retained for backward compatibility

### Function design
- Atomic functions first (detrend, window, filter, fit) — each independently useful
- Convenience wrappers second (analyze_diel_flux, analyze_annual_flux)
- Standard R inputs/outputs (zoo, tibble, numeric vectors) — no proprietary classes where avoidable
- Pipe-friendly signatures: data as first argument

### Code style
- tidyverse pipes (%>%) for data transformation
- purrr for iteration (map, map_dfr, safely) — no for loops
- roxygen2 documentation with @description, @details, @param, @return, @examples
- Comments explain WHY, not WHAT — rationale and methodological decisions
