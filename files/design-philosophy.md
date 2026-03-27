# temperheic: Design Philosophy

## Purpose hierarchy

1. **Foundational science.** Get the 1D ADE methods right — Luce (2013, 2017), van Kampen (2022) — with robust signal processing for real-world, non-ideal temperature data.

2. **Applied methodology.** Translate the academic literature into practical tools for restoration ecologists, consulting hydrologists, fisheries biologists, and engineers. These practitioners need reliable flux estimates, not a literature review.

3. **Accessibility.** Clear documentation, sensible defaults, real-data vignettes that mirror actual field workflows. The science is invisible to the end user.

This applied orientation distinguishes temperheic from VFLUX, 1DTempPro, and iFLOW, which were built by and for academic researchers.

## Atomic function design (Wickham philosophy)

Every exported function does one thing well, is independently useful, and composes cleanly.

- `detrend_temperature()` is valuable to someone who never calls `thObservedSeries()`
- `window_temperature()` might end up in an unrelated time series workflow
- Build atoms first, convenience wrappers second
- Never lock functionality inside a monolithic pipeline
- Accept and return standard R objects (zoo, tibble, numeric vectors)
- Let users choose: high-level wrapper for the 80% case, atomic pipeline for the 20%

## Signal processing module

1. **Atomic functions** — `detrend_temperature()`, `window_temperature()`, `filter_bandpass()` — each independently useful
2. **Convenience wrappers** — `analyze_diel_flux()` etc. — orchestrate common workflows

The existing `thObservedSeries()` is the core inverse model. Signal processing functions prepare data for it; wrappers call it repeatedly and organize output.
