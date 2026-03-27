# van Kampen et al. (2022) — Quick Reference

**Full citation:** van Kampen, R., Schneidewind, U., Anibas, C., Bertagnoli, A., Tonina, D., Vandersteen, G., Luce, C., Krause, S., & van Berkel, M. (2022). LPMLEn–A Frequency Domain Method to Estimate Vertical Streambed Fluxes and Sediment Thermal Properties in Semi-Infinite and Bounded Domains. *Water Resources Research*, 58, e2021WR030886.

**Full PDF:** `refs/van_kampen_et_al_2022.pdf`

## What this paper is for
Advanced multi-sensor, multi-frequency estimation method with formal uncertainty. Planned for Phase 8 of temperheic. More statistically rigorous than the η-method but more complex to implement.

## Key innovations over Luce 2013

1. **Multi-sensor (n sensors)** — uses all available depth sensors simultaneously for joint parameter estimation, not just pairs
2. **Two domain models** — semi-infinite (M_SI) and bounded (M_BD) with lower boundary condition
3. **Maximum likelihood estimation** — optimal handling of noise; formal confidence intervals
4. **Local polynomial spectral estimation** — reduces spectral leakage compared to standard FFT
5. **Transient compensation** — handles temperature drifts/trends in the data

## Mathematical framework

### Transfer function approach
Solution in Laplace domain:
Θ(z,s) = G_U(z,s,θ) · U_U(s) + G_L(z,s,θ) · U_L(s)

Parameters to estimate: θ = {D, V} (diffusivity and advection velocity)

### Semi-infinite domain (M_SI)
G_U(z,s,θ) = exp(λ₁(s) · (z − z_U))
G_L(z,s,θ) = 0

where λ₁ = V/(2D) − √(V²/(4D²) + s/D)

### Bounded domain (M_BD)
Uses both upper and lower boundary sensors. Transfer functions involve hyperbolic functions — see full PDF for expressions.

### Cost function
Maximum likelihood minimization over frequency-domain residuals, weighted by noise covariance. Solved numerically (not closed-form like Luce 2013).

## Relationship to Luce 2013 η-method

| Aspect | Luce 2013 (η) | van Kampen 2022 (LPMLEn) |
|--------|---------------|--------------------------|
| Sensors required | 2 (pair) | n ≥ 2 (joint estimation) |
| Domain model | Semi-infinite only | Semi-infinite + bounded |
| Estimation | Closed-form (explicit) | Numerical optimization (MLE) |
| Uncertainty | Error propagation (σ_η, σ_v*) | Cramér-Rao bound from Fisher info |
| Spectral method | OLS or FFT | Local polynomial (LP) |
| Transient handling | Requires detrending | Built-in compensation |
| Complexity | Simple, fast | More complex, slower |
| When to use | Standard diel analysis, quick estimates | High-quality uncertainty, multi-sensor, bounded domains |

## temperheic implementation status

Phase 8 — not yet started. Requires:
- Transfer function representations (semi-infinite + bounded)
- Local polynomial spectral estimation
- Maximum likelihood cost function + optimizer
- Multi-sensor joint estimation
- Model selection: semi-infinite vs. bounded

## Key equations to implement (read full PDF for details)
- Transfer functions: Eqs. 4-8
- Cost function: Eq. 12
- Fisher information / Cramér-Rao: Eqs. 15-17
- Model selection criterion: Eq. 18
