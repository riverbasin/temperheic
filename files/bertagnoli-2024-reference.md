# Bertagnoli et al. (2024) — Quick Reference

**Full citation:** Bertagnoli, A., Luce, C., van Kampen, R., Schneidewind, U., van Berkel, M., Tranmer, A. W., Vandersteen, G., Krause, S., & Tonina, D. (2024). iFLOW: A Framework and GUI to Quantify Effective Thermal Diffusivity and Advection in Permeable Materials From Temperature Time Series. *Water Resources Research*, 60, e2024WR037370.

**Full PDF:** `refs/bertagnoli_et_al_2024.pdf`

## What this paper is for
Describes iFLOW, a Python-based framework for the same problem domain as temperheic. Useful for: benchmarking temperheic against an independent implementation, understanding the modular workflow architecture, and accessing the 8 synthetic test cases in Table 1.

## iFLOW's modular architecture (Figure 2)
Three fundamental processing steps — temperheic maps to these:

1. **Mathematical model selection** — steady state, transient (semi-infinite or bounded)
2. **Signal processing** — FFT (amplitude, phase, SNR) or local polynomial
3. **Parameter estimation** — Bredehoeft (steady), analytical (Luce-type), MLEn (van Kampen-type)

The key design insight: separating these steps lets users inspect intermediate results and catch errors at each stage. temperheic's atomic function philosophy aligns with this.

## Benchmarking test cases (Table 1)
Eight synthetic boundary condition scenarios for validation. Each uses known D and V, generates synthetic temperature data, then recovers parameters.

| Case | Upper BC | Lower BC | Domain | Flow |
|------|----------|----------|--------|------|
| 1 | Single freq | Zero gradient | Semi-infinite | Downwelling |
| 2 | Single freq | Zero gradient | Semi-infinite | Upwelling |
| 3 | Multi freq | Zero gradient | Semi-infinite | Downwelling |
| 4 | Multi freq | Zero gradient | Semi-infinite | Upwelling |
| 5 | Single freq | Fixed temp | Bounded | Downwelling |
| 6 | Single freq | Fixed temp | Bounded | Upwelling |
| 7 | Multi freq | Fixed temp | Bounded | Downwelling |
| 8 | Multi freq | Fixed temp | Bounded | Upwelling |

These are planned for temperheic Phase 4b benchmarking.

## Key differences from temperheic

| Aspect | iFLOW (Python) | temperheic (R) |
|--------|---------------|----------------|
| Language | Python | R |
| Interface | GUI (PyQt) | Functions + planned Shiny app |
| Audience | Academic researchers | Applied practitioners |
| Methods | Bredehoeft, analytical, MLEn | Luce η-method (core), LPMLEn (planned) |
| Signal processing | FFT, local polynomial | OLS, FFT, detrending, windowing |
| Unique features | Interactive GUI, step-by-step wizard | Pipe-friendly API, tidyverse integration, convenience wrappers |

## What to read in the full PDF
- Figure 2: The modular workflow diagram — excellent conceptual reference
- Table 1: The 8 synthetic test cases with parameter values
- Section 3: Mathematical model descriptions (complements van Kampen 2022)
- Section 5: Performance results on synthetic data (benchmarking targets)
