# Power-law fits to the observed extremes: |δ| = A·V^−α, with V in Mpc³
These are least-squares fits of log|δ| against log V. The ± values are 1σ standard errors from only 3–5 points, so treat them as rough.
Script: fits.py. Results: fits.json. Chart: density_vs_scale_fits.png (make_fits_chart.py).

| fit | points | α | log10 A | scatter |
|-|-|-|-|-|
| Over-dense (chart line) | Coma, Shapley core, Saraswati, SGW | 0.67 ± 0.08 | 4.00 | 0.19 dex |
| Over-dense + Shapley whole | 5 | 0.70 ± 0.08 | 4.06 | 0.23 dex |
| Over-dense, no Coma | Shapley core, Saraswati, SGW, Shapley whole | 0.57 ± 0.16 | 3.33 | 0.21 dex |
| Superclusters only | Saraswati, SGW, Shapley whole | 0.13 ± 0.27 (unconstrained) | 0.75 | 0.12 dex |
| Voids (chart line) | Boötes, SDSS, Eridanus | 0.22 ± 0.10 | 0.95 | 0.12 dex |
| Voids + Local Hole | 4 | 0.19 ± 0.07 | 0.80 | 0.12 dex |
| Voids, saturating form, |δ|/(1−|δ|) = A·V^−α | 3 | 0.35 ± 0.23 | 2.06 | 0.28 dex |
| Saturating form + Local Hole | 4 | 0.32 ± 0.15 | 1.90 | 0.25 dex |

How fast standard cosmology's expected spread falls with volume (σ ∝ V^−s, linear top-hat):
- s = 0.42 at R = 100 Mpc
- s = 0.53 at R = 300 Mpc
- s = 0.59 at R = 1 Gpc
- s = 0.61 at R = 2 Gpc
- s = 0.62 at R = 4 Gpc
- For comparison, purely random scatter (Poisson / white noise) gives s = 0.5.
