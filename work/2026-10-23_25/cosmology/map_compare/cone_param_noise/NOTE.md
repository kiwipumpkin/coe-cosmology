# Cone-param noise: A/A_i with variance INSIDE H_obs

**Created:** 2026-09-18 06:33 PM CDT (UTC-0500)

Two families — both kept. Noise is **never** added by shuffling finished map pixels.

## Why ≠ pixel patchwork

| Wrong (patchwork) | This run |
|---|---|
| Build finished $A$, then mix/jitter pixels | Recompute $H_{\mathrm{obs}}$ with noisy $u$ or $d$ |
| Noise unrelated to CoE geometry | Velocity noise **along $\hat{r}$ from CoE** |
| Average already-noised sky pixels | Objects (or angles) evaluated natively |

## Family 1 — CONE-NATIVE (primary)

For each cone center $\hat{c}$ on the 72×35 grid from `cone75_maps.npz`:

1. Draw $N_{\mathrm{samp}}=64$ sightlines **uniformly in solid angle** within θ=75.0° of $\hat{c}$.
2. For each object $\hat{n}_i$: compute model $\mathbf{u}$, optionally add radial velocity noise and/or distance noise, form $H_{\mathrm{obs}}$ and $A_{\mathrm{raw}}$.
3. Cone $A$ = mean of $A_{\mathrm{raw}}$ over the $N_{\mathrm{samp}}$ objects.
4. Mean-renormalize the cone map → $A/A_i$.

Shared sightline draws: seed `42`. Noise seeds: vel=43,
dist=44, both-v=45, both-d=46.

Cards: `CARD_cone75_native_*.png`
Strip: `STRIP_conenative_baseline_vel_dist_both_invsq.png`

## Family 2 — ANGLE (spreadsheet-style companion)

One evaluation per cone-center direction $\hat{n}$ (the spreadsheet idea:
add random to velocity/distance **in the calculation** at that angle). Same
physics dials and σ choices. Independent seeds (100+ offsets vs cone-native).

Cards: `CARD_cone75_angle_*.png`
Strip: `STRIP_angle_baseline_vel_dist_both_invsq.png`

## Variants (each family)

| map | rule |
|-----|------|
| **baseline** | $d=d_0=\varepsilon D$, pure model $u$ |
| **vel-only** | $d=d_0$; $\mathbf{u}_r \leftarrow \mathbf{u}_r + \xi\sigma_u\hat{r}$ |
| **dist-only** | $d_i=d_0(1+\eta\sigma_{d,\mathrm{rel}})$, clip $>0.02\,d_0$; model $u$ |
| **BOTH** | vel + dist noise together on the same objects/angles |

### Rationalized σ

- **$\sigma_u = 0.05\,H\,d_0$** (primary; also 0.10 louder): fraction of shell Hubble
  velocity $V_H=H d_0$. Mild few-percent hit on $H_{\mathrm{obs}}$; order-of-magnitude
  akin to few×10²–10³ km/s if mapped to local-cluster scales (illustrative, not
  a cosmographic conversion). Noise is **radial from CoE** ($\hat{r}$).
- **$\sigma_{d,\mathrm{rel}} = 0.15$** (primary; also 0.25): ~10–20% depth/measurement
  mix is plausible for Migkas-like wide-$z$ samples within a large cone.

## Physics dials
- twoscale local-match, ρ(D)=Ω_m ρ_crit
- D/L=5.9, f_u=0.94, α=1.0, ε=0.1
- κ ≈ 0.230588; units H=G=D=1
- CoE $+x$ → Migkas $(\ell,b)=(280.0^\circ,-20.0^\circ)$

## Legacy sanity vs `cone75_maps.npz`

Legacy = fine full-sky $H_{\mathrm{obs}}$→$A$→θ=75° solid-angle cone average.

| comparison | max\|Δ\| |
|---|---|
| legacy vs **cone-native** baseline | 0.07338 |
| legacy vs cone75_ref | 0.00002 |
| cone-native baseline vs ref | 0.07338 |
| legacy vs **angle** baseline | 0.29685 |

Cone-native no-noise should be **close** to legacy/ref (same cone idea; Monte Carlo
vs pixel integral → small residual). Angle baseline is a different operator (point
sample) — larger Δ expected.

## Plotting
- Migkas Galactic: ℓ LEFT (`lon_disp=−ℓ_wrapped`)
- magma [0.85, 1.20], cards ~1400×900
- Compare grid: `GRID_conenative_vs_angle_invsq.png`

## A/A_i spans (invsq, mean-renormed)

### Cone-native
- **cn/baseline**: [0.8879, 1.0826]
- **cn/vel_s05**: [0.8926, 1.0837]
- **cn/vel_s10**: [0.8776, 1.0789]
- **cn/dist_s15**: [0.8882, 1.0828]
- **cn/dist_s25**: [0.8877, 1.0830]
- **cn/both_s05_s15**: [0.8941, 1.0860]

### Angle
- **ang/baseline**: [0.6521, 1.1657]
- **ang/vel_s05**: [0.5418, 1.2639]
- **ang/vel_s10**: [0.4532, 1.6586]
- **ang/dist_s15**: [0.6495, 1.1663]
- **ang/dist_s25**: [0.6422, 1.1660]
- **ang/both_s05_s15**: [0.5551, 1.2618]

## Files
- `run_cone_param_noise.py`, `cone_param_noise_maps.npz`, `manifest.json`, this NOTE


## Noise averaging (important)

In the **cone-native** family each cone mean averages $N=64$ independent object
noises, so measurement-error scatter on $A/A_i$ is suppressed by $\sim 1/\sqrt{N}$
relative to a single sightline. Maps stay close to the no-noise baseline unless
$\sigma$ is large — this is correct for “many objects per cone,” not a bug.

The **angle** family has one draw per cone-center, so the same $\sigma_u$ produces
a much larger pixel-to-pixel scatter (see compare grid). Use angle maps to see
the raw per-measurement effect; use cone-native for Migkas-like cone means.
