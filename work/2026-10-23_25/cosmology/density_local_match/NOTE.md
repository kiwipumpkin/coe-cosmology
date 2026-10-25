# Local-density match note

**Date:** 2026-09-18 (America/Chicago)

## What "local match" means

Chris asked whether the density **where the observer exists** matches the observed cosmic matter density.
Previously the binding sheet fixed **mean** enclosed density via κ=Ω_m/2=0.15 (so ρ_mean(<D)/ρ_crit=Ω_m),
but for a decreasing Phase A exponential at D/L~6–7 the **point density** ρ(D) was ~15× below Ω_m ρ_crit.

This run **inverts** that priority:

**ρ(D) = Ω_m ρ_crit = 0.3 ρ_crit** (by construction; verified 0.3000).

We still report ρ_mean(<D) and the κ implied by M_enc.

## κ vs Ω_m/2 tension (honest)

For any centrally decreasing profile with ρ(D)=Ω_m ρ_c:

ρ_mean(<D) > ρ(D) ⇒ κ = ρ_mean/(2 ρ_c) > Ω_m/2 = 0.15.

**Do not fake κ=0.15.** Recommended κ≈0.231 (>0.15).

## Recommended profile

| Quantity | Value |
|---|---|
| Profile | two-scale: ρ=ρ_core0 e^{-r/L}+ρ_u |
| D/L (core) | 5.9000 |
| L/D | 0.1695 |
| f_u | 0.9400 |
| ρ(D)/ρ_c | 0.3000 |
| ρ_mean(<D)/ρ_c | 0.4612 |
| κ_total | **0.2306** |
| κ_core / κ_u | 0.0896 / 0.1410 |
| AF mapping | `hratio` (H_obs/H; sheet-like) |
| R²_A / R²_B (balanced) | **0.6852 / 0.6855** |
| Amp, φ | +0.0700, 212.82° |
| R² sumopt (same profile) | 0.6694 / 0.7361 |
| R² without rotation | 0.4018 / 0.4334 |

Chris sheet ballpark ~0.74/0.73. Here under local-density: **~0.685/0.685** balanced (both ≥0.65);
sumopt ~0.67/0.74. Strict dual ≥0.70 not quite reached.

## Why not pure exponential alone?

- **invsq** A∝1/(H_obs/H)² under ρ(D)=Ω_m ρ_c: best min(R²)≈0.574 at D/L=6.70, κ=2.34 (κ inflated by local-match on a steep tail).
- **hratio** pure-exp: best ≈0.680/0.690 at D/L=1.50, κ=0.228 — competitive, but core is less CoE-concentrated than sheet D/L~5.5.
- **Two-scale** keeps sheet-like D/L~6, lifts ρ(D) with a floor, keeps κ~0.23 and R²~0.685/0.685.

Plummer-like scans did not beat two-scale under the density constraint.

## Migkas diagnostic

Sheet-style Amp×sin (not free A cos φ dipole; **not** a B-null claim). Digitised Fig.6 with error bars.
See `migkas_overlay_local_density.png`, `migkas_scatter_local_density.png`.

## Regenerable

```bash
/workspace/cosmology/.venv/bin/python /workspace/cosmology/density_local_match/run_density_local_match.py
```
