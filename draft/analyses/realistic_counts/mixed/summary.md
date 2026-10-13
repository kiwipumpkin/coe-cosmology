# Lobes built from a cluster + void mix

Chris: "They need to be mixed in some range of ratios." Every lobe holds both clusters and voids. The only question is how far its cluster:void ratio may differ from the cosmic mean.

## Bookkeeping (from `../../local_web/volumes_counts.md`)
- **Cosmic mean:** 6.6e8 clusters to 9.8e7 voids, which is **6.77 clusters per void** (r₀).
- **What one "cluster" and one "void" carry:**
  - A cluster is a ~2.3e14 M☉ knot plus its share of filaments and walls: +8.6e44 kg.
  - A void has R ≈ 29 Mpc and δ −0.7: −5.8e45 kg. Voids fill 78% of the volume.
  - At r₀ the two exactly balance, so the mean web is net zero.
- **What a cell's mix does:** clusters ×(1+a) and voids ×(1+b) give δ = 0.546(a − b) and a ratio r/r₀ = (1+a)/(1+b).
- **Main rule:** voids may not fill more than they already do (78%). So within a ratio deviation X the cell can reach −0.546·X/(1+X) ≤ δ ≤ 0.546·X.
  - In practice the void count stays fixed and the cluster count shifts.
  - A **1% ratio shift ≈ 0.55% density contrast.**
- **Solver:** box-constrained least squares (`mix_core.py`, `scan.py`) on a 0.5 Gpc grid. Same kernel and residual as before, net mass zero, no extra smoothing.
  - The counts at 1× use the smallest change that still reaches 1×, by adding L2 damping (`counts1x.py`).
- **Two cases:**
  - **Inside, mirror allowed:** 0.5–13.8 Gpc.
  - **Inside + outside, no mirror:** one mirror half inside, the other half outside at 14.8–30 Gpc.
  - **Correction to the last report:** that "staggered" field put one mirror half inside and the other half outside. Each half had both dense and empty parts. It was not "dense side inside, empty side outside."

## Result (`leftover_vs_ratio_deviation.png`)
|allowed ratio deviation|inside, mirror allowed|inside + outside, no mirror|
|-|-|-|
|±1%|185×|–|
|±3%|93×|91×|
|±5%|26×|54×|
|±6%|4.9×|–|
|**±6.5%**|**0.39× (reaches 1×)**|–|
|±10%|reached|26×|
|±20%|reached|4.8×|
|±50%|reached|1.5× (upper bound, solver not converged)|
|**±60%**|reached|**reached**|
|0.1×–10× of the mean|reached|reached|

**Minimum ratio deviation for 1×:**
- **Inside, mirror allowed: about ±6.2–6.5%** (6.0% gives 4.9×; 6.5% gives 0.39×). Lobes run from 6.44 to 7.12 clusters per void, against 6.77.
- **No mirror (inside + outside): somewhere between ±20% and ±60%.**
  - ±60% clearly works. ±20% gives 4.8×.
  - The 30–50% runs didn't fully converge (the 30% point at 6.7× is out of line for that reason), so the true threshold may be nearer 30–50%.

## Objects needed at 1×
Counts are clusters added to dense cells plus clusters removed from empty cells, with voids unchanged.

|case|allowed deviation|clusters to add/remove|objects already in the active volume|active volume (|δ| > 0.5%)|peak δ|
|-|-|-|-|-|-|
|inside, mirror allowed|unconstrained (0.1×–10×)|3.1e6|9.3e7|1,500 Gpc³|17%|
| |±20%|3.3e6|9.7e7|1,550 Gpc³|11%|
| |±10%|4.3e6|1.3e8|2,040 Gpc³|5.5%|
| |**±6.5% (threshold)**|**1.07e7**|3.2e8|**5,080 Gpc³**|3.5% (at the limit)|
|inside + outside, no mirror|unconstrained|7.6e6|1.8e8|2,860 Gpc³|42%|
| |**±60% (threshold)**|**9.0e6**|2.1e8|3,390 Gpc³|32%|

- **The change is small next to what's there:** 3–5% of the objects already in the active volume. At the threshold, ±1.07e7 changes against 3.2e8 present. For scale, the whole visible sphere has 6.6e8 clusters.
- **Swapping instead of adding:** if clusters and voids are swapped one for one (fixed total count), the inside threshold needs about 1.4e6 swaps. But that pushes void filling to 82%, slightly over the 78% limit.

## Does a tighter ratio limit need more lobes or wider ones? Wider.
- **Same eight lobes** at the same longitudes:
  - Dense (+) at ℓ ≈ 15°, 128°, 195° and 263°.
  - Empty (−) at ℓ ≈ 75°, 162°, 228° and 315°.
- **They grow** as the limit tightens. Active volume goes from 1,500 Gpc³ (no limit) to 2,040 (±10%) to 5,080 Gpc³ (±6.5%). That is 3.4× the volume, saturated at the limit.
  - The lobes thicken inward to about 9 Gpc and fill whole arcs of the ring.
  - The number of objects to change rises 3.4× (3.1e6 → 1.07e7).
- **Below about ±6%,** no amount of widening helps. Lobes far from the ring have too little lever on the CMB.

## Best mixed layout (`map_mixed_ratio.png`): inside, ±6.5%
|lobe ℓ|clusters per void|vs mean|mean δ|clusters added/removed|objects present|
|-|-|-|-|-|-|
|15°|7.01|+3.6%|+1.4%|+6.0e5|2.7e7|
|128°|7.12|+5.1%|+2.3%|+1.1e6|3.1e7|
|195°|7.12|+5.2%|+2.3%|+1.2e6|3.3e7|
|263°|7.12|+5.2%|+2.3%|+1.2e6|3.3e7|
|75°|6.48|−4.2%|−1.7%|−8.8e5|3.3e7|
|162°|6.44|−4.9%|−2.1%|−1.2e6|3.4e7|
|228°|6.44|−4.8%|−2.1%|−1.2e6|3.6e7|
|315°|6.48|−4.2%|−1.7%|−8.8e5|3.2e7|

- **By volume:** voids stay at 78% everywhere.
- **By mass:** the clusters' share of the mass moves from 28.6% to about 30.1% in dense lobes and about 27.2% in empty lobes.

## Can it be built?
- **Packing: yes, easily.** It's a ±4–5% shift in cluster counts with voids untouched.
- **Density: no, by observed standards.**
  - Inside, the lobes are coherent ±1.4–2.3% means (cells up to ±3.5%) over about 5,000 Gpc³, in a fixed m = 3–5 pattern 9–14 Gpc from us.
  - That is at or above the roughly 1% Gpc-scale homogeneity Chris quoted, and 20–200× the ΛCDM rms at those scales (see `../../inside_visible/smooth_field`).
  - The no-mirror version needs 30–60% ratio swings and δ up to ±32%. That is far outside anything observed.

## Caveats
- **No smoothing:** these fields are on a 0.5 Gpc grid. The earlier smooth-field run used 2 Gpc smoothing; smoothing would need slightly more.
- **Solver convergence:** the box solver didn't fully converge at small X. Results there are upper bounds (e.g. ±5%: 29× at 400 iterations, 26× at 1,500). So the thresholds could be slightly lower.
- **Threshold counts are not minimal:** the threshold solutions are undamped, so the count at the threshold (1.07e7) is an upper estimate.
- **Bookkeeping is the foam's:** cluster = knot + filament share, void = δ −0.7. With C&J-like voids (δ −0.4, filling 45–70%), voids could also be added, which slightly relaxes the bounds.
- **Missing physics (as before):** time-varying potentials, lensing and the ring–us potential difference are left out.
- **What was cut:**
  - The F = 1.28 variant (voids allowed to fill 100%).
  - Converged runs at 25–50% for the no-mirror case.
  - Damped minimal counts at the two thresholds.

## Files
- **Charts:** `leftover_vs_ratio_deviation.png`, `map_mixed_ratio.png`
- **Data:** `scan_results.json`, `counts_at_1x.json`, `lobes_mixed.json`, `sol_*.npz`, `best1x_*.npz`
- **Code:** `mix_core.py`, `scan.py`, `counts1x.py`, `report.py`
