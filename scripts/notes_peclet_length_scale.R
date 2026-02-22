
# Peclet Number Length Scale Discussion
# Date: 2026-02-22
# Context: temperheic development, dispersivity/diffusivity exploration
#
# The Peclet number Pe = vL/D is conceptually simple but practically slippery
# because the choice of length scale L determines what question you are asking.
#
# Three natural length scales appear in this work:
#
# 1. Grain diameter (d50) — pore-scale, used by Bons et al. (2013)
#    - Asks: is transport through individual pore channels dominated by
#      advection or diffusion?
#    - Critical Pe_c ~ 3-40 marks transition from diffusion-dominated to
#      mechanical dispersion regime at the pore scale.
#    - Material property, independent of measurement setup.
#
# 2. Thermal decay distance (x_d) — signal-scale, emerges from Vogt (2012)
#    - The e-folding distance for amplitude damping: x_d = f(omega, v_t, kappa_e)
#    - Not a material property — depends on the frequency of the signal being
#      tracked (diel vs annual) and the flow/diffusion regime.
#    - Diel x_d is O(1 m) in streambeds; annual x_d is O(100-1000 m) in
#      floodplains.
#    - Used implicitly in Luce (2013) nondimensionalization: x = z/z_d.
#
# 3. Sensor/well spacing — measurement-scale
#    - Asks: does advection or diffusion dominate transport over the distance
#      between my observation points?
#    - Practical but arbitrary — depends on field setup.
#
# Key insight: the ratio kappa_e / kappa_cond is a dimensionless diagnostic
# that avoids the length scale choice entirely. It directly quantifies how
# much of the effective thermal transport is due to mechanical dispersion
# vs conduction:
#   kappa_e / kappa_cond ~ 1  -> conduction-dominated (low Pe regime)
#   kappa_e / kappa_cond >> 1 -> dispersion-dominated (high Pe regime)
#
# This ratio is recoverable from the inverse model (Luce 2013) without
# committing to a length scale for Pe.
#
# The multi-frequency approach (diel + annual eta-consistency, Luce 2017)
# naturally probes different effective length scales because diel and annual
# signals have different x_d. Comparing eta across frequencies implicitly
# tests whether the Peclet regime is consistent across scales — the Bons
# question without needing to pick a single L.
#
# The "diffusive" and "dispersive" Peclet decomposition explored in the
# original pecletNumbers.R script was probing this same idea: partitioning
# effective transport into conductive and dispersive contributions. While
# not standard nomenclature in the literature, it captures a real physical
# distinction. The kappa_e/kappa_cond ratio and multi-frequency eta
# consistency are more rigorous ways to get at the same question.
#
# References:
#   Bons et al. (2013) WRR — unified D(Pe) expression, Pe regimes
#   Vogt et al. (2012) HESS — x_d / v_d formulation
#   Luce et al. (2013) WRR — eta parameter, nondimensional z/z_d
#   Luce et al. (2017) — multi-frequency eta consistency diagnostic

