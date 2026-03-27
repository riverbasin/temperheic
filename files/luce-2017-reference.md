# Luce et al. (2017) — Quick Reference

**Full citation:** Luce, C. H., Tonina, D., Applebee, R., & DeWeese, T. (2017). Was That Assumption Necessary? Reconsidering Boundary Conditions for Analytical Solutions to Estimate Streambed Fluxes. *Water Resources Research*, 53, 9771-9790.

**Full PDF:** `refs/luce_et_al_2017.pdf`

## What this paper is for
Proves that the Luce 2013 η-method works even when the surface temperature signal isn't a clean sinusoid. Eliminates two common objections to the analytical approach. Provides the multi-frequency consistency diagnostic.

## Key findings

### 1. Boundary condition independence
The η-method does NOT require:
- A perfect sinusoidal surface boundary condition
- Zero mean temperature gradient with depth

The parameter estimates are valid for any arbitrary Fourier-decomposable surface signal. This means real-world noisy, non-sinusoidal stream temperature records work fine.

### 2. Frequency independence
For any frequency ω in the signal:

v* = (1 − η_ω²) / √(2(1 + η_ω²))

κ_e = ω · Δz² / (Δφ_ω² · (1/η_ω + η_ω))

Each frequency should independently give the same v* and κ_e estimates. If they don't, something is wrong (heterogeneity, violated assumptions, bad data).

### 3. Multi-frequency consistency diagnostic
Compare η_ω across multiple frequencies (diel, semi-diurnal, harmonics). Consistent η values validate the homogeneous 1D assumption. Divergent values indicate:
- Vertical heterogeneity in thermal properties
- Lateral flow contributions
- Sensor problems
- Non-1D geometry

## Practical implications for temperheic
- `fit_ols()` and `fit_fft()` can target any frequency, not just the fundamental diel
- Multi-frequency η comparison is a powerful diagnostic (planned for Phase 2b Example 6)
- Users don't need to pre-filter their data to isolate a clean sinusoid
- Annual and diel analyses of the same record should yield consistent flux estimates (within measurement uncertainty) if assumptions hold

## temperheic implementation mapping

| Paper concept | temperheic status | Notes |
|---------------|------------------|-------|
| Single-frequency η at arbitrary ω | Implemented in `fit_ols()`, `fit_fft()` | `period` parameter |
| Multi-frequency consistency check | Not yet implemented | Planned: Phase 2b Example 6 |
| Boundary condition generality | Implicitly supported | OLS fitting handles arbitrary signals |
