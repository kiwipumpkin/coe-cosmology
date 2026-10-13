# Density distribution vs region size: summary

Charts:
- density_hist_by_size.png: how common each δ is, for 10, 20, 40, 80 and 150 Mpc spheres.
- density_vs_size_heatmap.png: the same axes as the v2 chart, with typical regions as a heat map under the extremes.
- object_counts_vs_size.png: number of clusters, superclusters and voids per Gpc³ vs radius.

Numbers are in distribution_stats.json, cic_2mpp.json and object_counts.json.

## 1. Shape and width vs size (2M++, mass δ)
| R (Mpc) | std | ΛCDM linear σ (same smoothing) | 1st–99th pct | independent spheres |
|---|---|---|---|---|
| 10 | 0.54 | 0.65 | −0.62 … +2.1 | ~9,800 |
| 20 | 0.40 | 0.47 | … +1.4 | ~1,230 |
| 40 | 0.24 | 0.27 | −0.37 … +0.68 | ~150 |
| 80 | 0.10 | 0.13 | −0.17 … +0.25 | ~19 |
| 150 | 0.02 | 0.06 | ~±0.03 about +0.04 | ~1–3 |

- **Shape.** At small sizes the distribution is lopsided: under-dense regions stop at −1, and over-dense regions have a long tail (skew 2.6 at 10 Mpc). By 80 Mpc it is close to a symmetric bell curve (skew ≈ 0).
- **Width.** The 2M++ spread falls as std ≈ 8.6 · R^−1.00±0.09 (R in Mpc) over 20–80 Mpc.
- **Versus ΛCDM.** At 20–80 Mpc the 2M++ spread is about 0.82 × the ΛCDM value, which agrees within bias and smoothing uncertainty. At ≥110 Mpc 2M++ falls below ΛCDM, because there are too few independent regions and the reconstruction shrinks toward the mean.
- **Extremes.** The published extreme objects sit in the far tails of these histograms, at roughly 3–8 std (several lie beyond the densest sphere in the local 2M++ volume). That is where the extremes on the v2 chart belong: they are the rare tail, not the typical region.

## 2. How many objects of each size (per Gpc³ per factor-10 in radius)
- Clusters (redMaPPer, λ ≥ 20): peak ~2.4×10⁴ at R200 ≈ 1.8 Mpc, then a steep fall to ~20 at 4.5 Mpc. This follows the Tinker+2008 model.
- Superclusters (Liivamägi+2012): broad peak ~5×10³ at R ≈ 9 Mpc, falling to ~70 by 28 Mpc.
- Voids (BOSS, Mao+2017): peak ~230 at R ≈ 90 Mpc, falling to ~3 at 350 Mpc. The largest void has R ≈ 650 Mpc (453 h⁻¹).
  - The model (Jennings+2013) counts matter-defined voids, which are much smaller, so model and data do not overlap. That reflects the different definitions, not a contradiction.
- Nothing in any catalog reaches 1 Gpc.

## 3. Rarity of the required field (peak |δ| = 0.068 / 0.055 / 0.059 at 1 / 2 / 4 Gpc)
Two-sided Gaussian probability for one region. "Expected" is how many such regions the 14.3 Gpc visible sphere would hold.

| R | width used | σ | peak in σ | P (one region) | expected in visible sphere |
|---|---|---|---|---|---|
| 1 Gpc | 2M++ trend (generous) | 0.0085 (0.0065–0.011) | 8.0 (6.0–10.5) | 1e−15 (1e−9 to 1e−26) | 4e−12 (2,900 regions) |
| 1 Gpc | ΛCDM | 0.0028 | 24.5 | 1e−132 | ~0 |
| 2 Gpc | 2M++ trend | 0.0043 | 12.9 (9–18) | 3e−38 | 1e−35 |
| 2 Gpc | ΛCDM | 0.00079 | 69 | ~1e−1050 | 0 |
| 4 Gpc | 2M++ trend | 0.0021 | 28 (19–41) | 1e−169 | 0 |
| 4 Gpc | ΛCDM | 0.00022 | 267 | ~1e−15500 | 0 |

- The required *rms* (0.0105 / 0.0124 / 0.0141) is 1.2× / 2.9× / 6.6× the generous 2M++ trend, and 3.8× / 16× / 64× ΛCDM.
- Volume fraction with |δ| > 0.055 in 2M++: 0.82 (10 Mpc), 0.85 (20), 0.80 (40), 0.62 (80), 0.17 (150; this is mostly the local mean offset of +0.04). The fraction is already dropping fast by 150 Mpc.
- **Bottom line.** A ±5–7% field on Gpc scales is about 8σ at 1 Gpc even with the most generous width the data allow, and 25–270σ in ΛCDM. The model therefore cannot be a "typical" region drawn from the observed distribution. It would need a new, much larger source of large-scale variation; required rms ÷ the ΛCDM width is 4–64×.

## Caveats
- The extrapolation from 80 Mpc to 1–4 Gpc covers a factor of 12–50 in radius.
- The 2M++ slope (−1.0) is shallower than ΛCDM's (−1.25 to −1.6 at these sizes). That makes the extrapolation generous to the required field.
- Gaussian tails are assumed at Gpc scales (appropriate when σ ≪ 1).
- See sources.md for derived and approximate flags.
