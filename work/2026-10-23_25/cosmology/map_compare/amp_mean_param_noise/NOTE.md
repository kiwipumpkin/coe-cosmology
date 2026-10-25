# Amp×sin family mean base + experimental v/d noise

**Created:** 2026-09-18 06:39 PM CDT

**Base:** exact stored mean of 8 Amp×sin cone maps (`amp_family_mean.npz` /
`CARD_amp_family_MEAN.png`), span ~[0.916, 1.028] before any re-renorm (mean 1).

**Noise (experimental):** at each θ=75° cone center, recompute model \(H_\mathrm{obs}\)
with radial velocity noise and/or distance noise; form
`factor = A_invsq(noisy) / A_invsq(true)`; set `A = mean_renorm(A_base * factor)`.

This keeps your base morphology and adds only the measurement scatter.

- σ_u = 0.05 H d_0 (radial from CoE)
- σ_d,rel = 0.15

Spans:
- baseline: [0.9162, 1.0277]
- vel: [0.7402, 1.2871]
- dist: [0.9129, 1.0311]
- both: [0.6545, 1.5332]
