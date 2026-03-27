# Thermal Hydrology Domain Conventions (Always Active)

## Physical parameters — always verify units

| Parameter | Symbol | Typical range | Units |
|-----------|--------|---------------|-------|
| Hydraulic conductivity | K | 1e-7 to 1e-2 | m/s |
| Thermal diffusivity | κ_e | 1e-7 to 1e-6 | m²/s |
| Thermal velocity | v_t | -1e-5 to 1e-5 | m/s |
| Darcy flux | q | -1e-5 to 1e-5 | m/s |
| Bulk thermal conductivity | λ_m | 0.5 to 3.0 | W/(m·°C) |
| Bulk heat capacity | ρc | 1.5e6 to 4.0e6 | J/(m³·°C) |
| Porosity | n | 0.2 to 0.5 | dimensionless |
| Dispersivity | ψ | 0.001 to 0.1 | m |

## Key relationships

- κ_e = λ_m / (ρ_m · c_m) — effective thermal diffusivity
- v_t = q · ρ_w · c_w / (ρ_m · c_m) — thermal front velocity from Darcy flux
- η = ln(A₂/A₁) / (φ₂ - φ₁) — ratio of amplitude damping to phase shift
- v* = (1 - η²) / √(2(1 + η²)) — dimensionless velocity from η
- η < 1 → downwelling; η = 1 → no flow; η > 1 → upwelling

## Signal characteristics

- Diel signal: period = 86400 s; penetrates ~0.5-2 m in streambeds
- Annual signal: period = 31557600 s; penetrates ~10-100 m in floodplains
- Semi-diurnal (12h) and higher harmonics exist but are weaker
- Amplitude attenuates exponentially with depth; phase lags linearly (in pure diffusion)

## Meacham Creek study site context

- Floodplain restoration site on Meacham Creek, NE Oregon
- Stream + 17 wells, 2011-2016 monitoring
- Key pairs: Stream–Well 20 (shallow, diel), Stream–Well 21–Well 10 (transect, annual)
- Air temperature proxy from UP211 SNOTEL station (1.5 km downstream)
- Surface water record reconstructed via regression against air temperature
