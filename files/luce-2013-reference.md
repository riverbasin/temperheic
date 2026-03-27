# Luce et al. (2013) — Quick Reference

**Full citation:** Luce, C. H., Tonina, D., Gariglio, F., & Applebee, R. (2013). Solutions for the diurnally forced advection-diffusion equation to estimate bulk fluid velocity and diffusivity in streambeds from temperature time series. *Water Resources Research*, 49, 488-506.

**Full PDF:** `refs/luce_et_al_2013.pdf` (read for exact derivations, figures, sensitivity contour plots)

## What this paper is for
The foundational method implemented in temperheic. Provides closed-form analytical solutions for estimating groundwater flux and thermal diffusivity from temperature time series at two depths. The η parameter is the core innovation.

## Governing equation (1D advection-diffusion for heat)

∂T/∂t = κ_e · ∂²T/∂z² − v_t · ∂T/∂z

- κ_e = λ_m / (ρ_m · c_m) — effective thermal diffusivity [m²/s]
- v_t = q · ρ_w · c_w / (ρ_m · c_m) — thermal front velocity [m/s]

## The η parameter (Eq. 5 in paper)

η = ln(A₂/A₁) / (φ₂ − φ₁) = α/β

- A₁, A₂: amplitudes at depths z₁, z₂
- φ₁, φ₂: phases at depths z₁, z₂ (radians)
- η < 1 → downwelling (v* > 0)
- η = 1 → no flow (v* = 0)
- η > 1 → upwelling (v* < 0)

## Dimensionless velocity from η (Eq. 6)

v* = (1 − η²) / √(2(1 + η²))

This is the key explicit solution — no iteration needed.

## Inverse: η from v* (Eq. 7)

η = (√(v*⁴ + 4 + v*²) − √2 · v*) / √(v*⁴ + 4 − v*²)

Implemented in `thSeries()` and `thObservedSeries()`.

## Thermal diffusivity from known depths (Eq. 8)

κ_e = ω · Δz² / (ln²(A_r) + Δφ²)

Equivalent form: κ_e = ω · Δz² / (Δφ² · (1/η + η))

## Thermal velocity (dimensional) (Eq. 9)

v_t = (ω · Δz / Δφ) · (1 − η²) / (1 + η²)

## Nondimensional scaling

- z_d = √(2κ_e/ω) — damping depth
- v_d = √(2ωκ_e) — diffusive phase velocity
- v* = v_t / v_d — dimensionless velocity

## Error propagation — UQ equations (Gariglio et al. 2013)

### Uncertainty in η (Eq. 13 in Gariglio)
σ_η = √[ σ_A² / (φ₂−φ₁)² · (1/A₂² + 1/A₁²) + 2·σ_φ²·ln²(A₂/A₁) / (φ₂−φ₁)⁴ ]

### Uncertainty in v* (Eq. 14 in Gariglio)
Complex expression — see full PDF or `Luce_et_al_2013_Streambed_Thermal_Analysis.md` in docs/

### Uncertainty in v_t from κ_e uncertainty
σ_{v_t,κ_e} = v* · √(ω / (2κ_e)) · σ_{κ_e}

## Validation approaches
1. κ_e consistency across seasons (should be stable despite changing T and flow)
2. Multi-frequency η consistency (Luce 2017 diagnostic — each ω should give same v*)
3. Streambed elevation tracking (Δz estimation when κ_e is known)

## temperheic implementation mapping

| Paper concept | temperheic function | Notes |
|---------------|-------------------|-------|
| OLS amplitude/phase extraction | `fit_ols()` | Uses linearized cosine: A·cos(ωt+φ) = α_s·sin + α_c·cos |
| η calculation | Internal to `thObservedSeries()` | |
| v* from η | `thObservedSeries()` | |
| κ_e from (η, Δz, ω) | `thObservedSeries()` | |
| Forward T(z,t) | `thSeries()` | Given v*, κ_e, boundary conditions |
| Round-trip validation | `test_round_trip()` | generate → fit → recover |
| UQ: σ_η, σ_v* | Phase 1 — not yet implemented | |
