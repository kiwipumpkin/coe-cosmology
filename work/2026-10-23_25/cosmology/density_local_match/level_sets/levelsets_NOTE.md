# Level-set guide: which dial moves what

Baseline (balanced recommended profile, held when not scanning that axis):

- two-scale, **D/L = 5.9**, **f_u = 0.94**, **κ ≈ 0.231**, AF = `hratio`
- Amp = **+0.07**, φ ≈ **213°** (3.714 rad), **α = 1**, **d/D = 0.1**
- Constraint: **ρ(D) = Ω_m ρ_crit = 0.3 ρ_crit** on every profile cell
- Amp_A, Amp_B independent with **shared φ**, re-fit per cell on profile panels
  (best of indep-amps vs shared-opposite under max-min R²)

## One-sentence dial map

| Dial | What it primarily changes |
|------|---------------------------|
| **D/L** | Core compactness of the exponential core → tidal amplitude of H_obs and **κ** |
| **f_u** | Uniform outer mass fraction → boosts κ with weaker core tide; trades shape |
| **α** | Gravity / peculiar-velocity **strength** (overall anisotropy scale) |
| **d/D (eps)** | Observer offset → **multipole mix** (dipole vs quadrupole in AF/AN) |
| **φ** | Sky **axis aim** of the Amp×sin template |
| **Amp_A / Amp_B** | Template **depth** on A vs B (phase shared) |

κ is reported as white dashed contours wherever the density profile changes
(D/L–f_u and D/L–α). κ is fixed on α–eps and Amp–φ panels (same baseline profile).

## Panel peaks (this run)

- **D/L vs f_u**: best min(R²) ≈ 0.692 at D/L=5.69, f_u=0.980, κ=0.174
  (R²_A=0.692, R²_B=0.747, Amp_A=+0.020, Amp_B=-0.150, φ=217.0°)
- **α vs d/D**: best min(R²) ≈ 0.710 at α=0.62, d/D=0.220
  (R²_A=0.710, R²_B=0.743)
- **D/L vs α**: best min(R²) ≈ 0.692 at D/L=5.69, α=0.50, κ=0.221

Baseline verification (ρ(D)/ρ_crit): 0.3000 (target 0.3).
Baseline re-fit R²_A=0.6854, R²_B=0.7467
(Amp_A=-0.070, Amp_B=+0.110, φ=34.5°).

## Files

- `levelsets_DL_vs_fu.png` (+ pdf)
- `levelsets_alpha_vs_eps.png` (+ pdf)
- `levelsets_Amp_vs_phi.png` (+ pdf)
- `levelsets_DL_vs_alpha.png` (+ pdf)
- `levelsets_summary.png` (+ pdf) — combined + caption strip
- `run_levelsets.py` — regenerable

Digitised Migkas: `/workspace/cosmology/migkas/fig6_digitized.csv` (bars not drawn on contours).
