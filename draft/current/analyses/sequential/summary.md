# Sequential, non-overlapping placement — Scenarios A and B

Every object is a sphere placed **one at a time**. No two spheres overlap: every pair is at least R_i + R_j + 0.02 Gpc apart, checked by brute force on every final layout (smallest gap found is 0.020 Gpc).
- **Voids:** δ = −0.70.
- **Clusters:** δ = +100, so they are compact, about 0.5% of the volume. Clusters carry exactly the mass the voids removed, so net mass is zero.
- **Scoring:** same tidal pattern, residual and 1e-5 floor as before. Start is 237×; the four clumps reach 89×.
- **Region:** 14.3–42.9 Gpc from us, with the CoE interior (<43.9 Gpc from the CoE) left empty.
- **Coarse-graining:** void radius is 1.5 Gpc in the outer fills and in B, and 1.5–3.75 Gpc in the near fit. Real voids are about 0.03 Gpc, so these are 50–100× larger.
- **Packing limit:** hard spheres can't fill 78% of space. The fills jam at about 30–34% void volume (lattice 20%).

## Scenario A — fix the near region, freeze it, then fill the outer shell
**Step 1, near region (14.3–30 Gpc):** 24 void + cluster pairs were placed one at a time, each pair picked to lower the leftover the most. Across 4 seeds the leftover went from 237× to between **27× and 45×**; the best is **27.4×**, well below the four clumps' 89×.
- Every mode is reduced, not just the dipole: m=1 went 18× → 1.5×, and m=3–5 went from about 200× to about 1–3×.
- The objects are supercluster/supervoid sized, comparable to the four clumps (voids 1.5–3.75 Gpc, clusters up to about 7e51 kg).

**Step 2, near region frozen at 27.4×; outer shell (30–42.9 Gpc) filled one unit at a time.** A unit is a void, then 6 clusters on its rim.

|outer fill|layouts|median leftover|range (10–90%)|helped / hurt|
|-|-|-|-|-|
|empty (control)|1|27.4×|—|—|
|lattice, ordered (~2,500 units)|30|27.4×|27.4–27.4|0 / 0 (change < 0.01×)|
|foam-like, grown next to neighbours (~3,900 units)|30|28.7×|28.1–30.0|0 / 30|
|random packing (~3,900 units)|30|29.2×|28.4–30.5|0 / 30|
|48 big void/cluster pairs scattered at random|200|152.8×|98–236|0 / 200|

**Answer:**
- Once the near region is set, a realistic, locally balanced outer web hardly matters: about +1 to +2× on 27×.
- An ordered lattice changes nothing.
- What does matter is big unbalanced structures (supercluster-scale pairs) in the outer shell. Those undo most of the fix.

## Scenario B — grow from empty, one object at a time (no preset fix)
Void radius is 1.5 Gpc. Growth runs until the region jams, at about 5,500–6,500 voids and 33,000–39,000 clusters (alternate: one heavy cluster per void).

|growth rule|seeds|lowest point while growing|when full|
|-|-|-|-|
|foam: voids next to voids, clusters at junctions|3|221–229×|438–624×|
|alternate: void, then one heavy cluster across a neighbour|3|214–237×|772–912×|
|CoE-biased: voids away from the CoE, clusters facing it|3|225–230×|671–849×|
|random: same objects dropped anywhere (control)|3|213–217×|1,045–1,621×|
|CMB-aware: each spot chosen to lower the leftover|2|36–39×|169–203×|

- Structure-only rules (foam, alternate, CoE-biased) never drift far below 237×. As the shell fills they all hurt.
- Structured growth is still clearly better than random scatter (foam about 440–620× against random about 1,000–1,600×).
- Only CMB-aware growth gets below the four clumps. It drops past 89× after a few hundred voids and reaches about 36×. It rises again at the very end, when jammed leftover spots force bad placements.
- The coarse 1.5 Gpc voids make the totals far too large. Earlier foam runs showed the web's own pattern shrinks roughly in proportion to void size, so real 0.03 Gpc voids would be about 50× smaller.

## Charts
- `a_local_fit.png`: near-region leftover vs pairs placed
- `a_local_layout.png`: the frozen near arrangement, to scale (nothing overlaps)
- `a_outer_effect.png`: outer-fill configs vs the frozen 27.4×
- `a_helped_hurt.png`: helped/hurt vs the frozen value
- `b_growth.png`, `b_growth_first300.png`: leftover vs voids placed
- `b_final_bars.png`: lowest point and final leftover for each rule
- `b_example_slices.png`: to-scale slices of foam vs random growth

## Cuts and limits
- Dropped: the outer fill with half-size voids (0.75 Gpc). It did not finish one layout in 40 minutes.
- Dropped: the B runs at 1.0 Gpc voids.
- CMB-aware growth ran 2 seeds; the other rules ran 3.
- Exact pair overlap checks were done on every layout.
