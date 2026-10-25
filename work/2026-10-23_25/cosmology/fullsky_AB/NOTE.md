# Full-sky A(n̂) and B(n̂) — first principles + sheet sine diagnostic

**Date:** 2026-09-18 01:18 PM CDT

## What A and B mean here (PRIMARY — first principles)

These maps are **physics-only** angular series from Chris’s ST+Gauss two-scale geometry
(binding baseline: `twoscale_fu0.94`, D/L=5.9, f_u=0.94, κ≈0.231, ρ(D)=0.3 ρ_crit).

| Symbol | Definition | Mean-renorm |
|---|---|---|
| **A_hratio** (primary AF) | A_raw = H_obs/H | A / ⟨A⟩ |
| **A_invsq** (secondary) | A_raw = 1/(H_obs/H)² (H_obs>0) | A / ⟨A⟩ |
| **B** | AN = m(r_src)/M_enc | B / ⟨B⟩ |
| **frac** | H_obs/H − 1 | (not renormed) |

**Honest:** B here is the **mass-proxy AN** used in `density_local_match`, **not** the Migkas
L_X–T slope parameter B from the observational fit. Do not conflate them.

Primary Mollweides (`mollweide_AB_*_frame.png`, `mollweide_sky_*.png`) use **raw physics
H_obs and m(r) only** — Amp×sin is **not** mixed into those maps.

## Sheet sine rotation — DIAGNOSTIC LAYER (retained per Chris)

Chris wants the unexplained sheet Amp×sin phase rotation **kept** as a labeled
diagnostic. It is **not** derived from Gauss / ST and is **not** first-principles gravity.

Sheet convention (applied as a sky function of Galactic longitude ℓ of each pixel):

```
BA = Amp_A · sin(2π ℓ/360 + φ)
BB = Amp_B · sin(2π ℓ/360 + φ)   # Amp_B = −Amp_A (opposite)
Model A = (1 − BA) · AF
Model B = (1 − BB) · AN
```

Two diagnostic variants (sky-frame AF=A_hratio, AN=B):

1. **sheet defaults:** Amp_A=−0.06, Amp_B=+0.06, φ=0.35 rad
2. **local-density balanced fit:** Amp_A=+0.07, Amp_B=−0.07, φ≈212.82°

Outputs clearly labeled:
`mollweide_sine_*_DIAGNOSTIC.png` — **"sheet sine rotation diagnostic — not derived from Gauss"**.

Do **not** present sine-rotated maps as first-principles gravity. They are an
axis/modulation diagnostic Chris wants retained for Migkas Fig.6-style comparisons.

## Quadrupole-led expectation (physics frac)

At d/D = ε = 0.1, the tidal limit of ST+Gauss gives Ψ/H ≈ ακ(3μ²−1) = 2ακ P₂(μ):
**pure quadrupole, vanishing dipole at leading order.** Finite ε generates a small
dipole. Multipole moments on `frac` (model frame):

- |A_dipole| = 4.8390e-03
- Q_P2 = 1.6012e-01
- |C1| ≪ |C2|: |A| / |Q_P2| = 3.0222e-02

## Frames and rotation

- **Model frame:** CoE at origin, observer at +D x̂, preferred axis D̂ = +x.
  Lon measured from +x in the xy-plane; lat = elevation.
- **Sky frame (display):** rotate so +x → Migkas-ish preferred direction
  (ℓ, b) ≈ (280.0°, -20.0°).
  Rotation R maps e_x → n_gal(ℓ,b); sky value at n_sky = model value at Rᵀ n_sky.
  Documented in `metadata.json` (`rotation_matrix`, `migkas_lb_deg`).
- Sine diagnostic uses **Galactic ℓ of each sky-frame pixel** (lon on the sky grid).

## Baseline parameters

- profile: twoscale, D/L=5.9, f_u=0.94
- κ = 0.230588, κ_core = 0.089588
- α = 1.0, d/D = 0.1, Ω_m = 0.3
- grid: n_lon=360, n_lat=181

## Next step

θ = 75° cone average (directional averages of A/B over Migkas-style cones)
for direct comparison to Fig.6-style longitude series. Cone averages can be run
on physics-only maps and, separately, on the labeled sine-diagnostic maps.

## Regenerable

```bash
/workspace/cosmology/.venv/bin/python /workspace/cosmology/fullsky_AB/run_fullsky_AB.py
```
