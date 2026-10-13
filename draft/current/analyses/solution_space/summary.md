# Fitting masses and voids on purpose: the solution space

The question: if masses and voids are placed deliberately to cancel the leftover CMB pattern (237× the 1e-5 floor), how many objects does it take, how many different arrangements work, and what do the working arrangements have in common?

## Setup
- **Target and score:** the same leftover R (237.3×), the same tidal kernel (`field3` / `web_core.pattern`) and the same floor as before. Leftover = rms(R − fitted pattern) / floor.
- **Region:** 14.3–43 Gpc from us in any direction, with no core inside our sphere and none inside the CoE interior (< 43.9 Gpc from the CoE).
  - The old four clumps at ℓ 162°/228° sat inside the CoE interior, so they are **not** allowed here.
- **Object types:**
  - **Mass:** a compact cluster (δ = +100 core), M from 1e49 to 2e52 kg.
  - **Deep void:** δ = −0.70, no compensation, radius 0.5–4 Gpc. Mass is set by radius and local density, up to ~1.6e52 kg.
  - **Shallow void:** the Clampitt & Jain profile (mean −0.35, compensated out to 3 Rv), radius 0.5–4 Gpc.
- **Net mass:** any net mass (masses minus deep-void deficits) is balanced by the same smooth compensation the four-clump baseline used. Shallow voids are already net zero.
  - Check: the old four clumps give 88.8× (compensated) in this code, matching the earlier value.
- **Overlap:** cores never overlap (gap ≥ 0.02 Gpc). Every reported solution passes this check, and 100% of starts were feasible.
- **Fit:** each object has 4 free numbers (distance, longitude, latitude, size or mass), all fitted together. Least squares was run from random starts.
- **Classes:** N = 1–8 objects. That gives masses only, deep voids only, shallow voids only, and every split of masses + deep voids and of masses + shallow voids, for **80 classes** in total.
- **Starts (cut for time):** 100 random starts for 13 classes (the N = 8 ones and masses-only N = 7, done first). The other 67 had 40 starts each. The first launch ran ~4× slower than estimated, so the rest were relaunched at 40.
  - "Distinct solutions" are counted using the **first 40 starts in every class**, so the counts are comparable.
- **Distinct solution:** after removing objects that do nothing (own pattern < 1× floor), two solutions are the same if every object matches within 2 Gpc and 0.25 dex in mass (0.5 Gpc in void radius).
  - Exact symmetries are merged: flipping any object above or below the ring plane, and the mirror about the ℓ = 15°/195° axis (the residual and the CoE direction are both symmetric about it).

## 1. How fast the leftover drops (`a_best_vs_N.png`)
Best leftover found (× floor). Mixed classes show the best split.

|N|masses only|masses + deep voids|masses + shallow voids|deep voids only|shallow voids only|
|-|-|-|-|-|-|
|1|209|209|209|210|223|
|2|164|160 (1+1)|164|164|217|
|3|121|92 (2+1)|121|117|217|
|4|**36**|75 (3+1)|121|109|217|
|5|29|**22** (4+1)|34|101|188|
|6|23|16 (4+2 / 2+4)|29|95|188|
|7|17|**7.4** (4+3)|23|90|188|
|8|13|**4.5** (4+4)|17|85|143|

- **Masses are required.** Deep voids alone never beat the four clumps' 89× (85× at N = 8). Most of them sit at the 4 Gpc size cap, so they can't get strong enough.
- **Shallow voids do almost nothing.** Alone, they stay at 143–223×. Added to masses, the fit leaves them inert: for example, 7 masses + 1 shallow void gives exactly the 7-mass result.
- **4 masses is the key step.** It drops 121× → 36×. This beats the 89× four clumps (they were limited to the ring plane and a shared distance) and stays outside the CoE interior.
- **Adding deep voids to ≥ 4 masses** gives the best results: 7.4× with 4 masses + 3 voids, and 4.5× with 4 + 4.
  - In the best fits the total void deficit nearly equals the total cluster mass (4.4e52 vs −4.4e52 kg at 4 + 3; 6.2e52 vs −6.3e52 at 4 + 4).
- **Nothing reached the floor.** The best is 4.5×, and what's left is in m = 9–12.

## 2. How many different arrangements work (`b_distinct_below27.png`, `b_distinct_below50.png`)
Distinct solutions in the first 40 starts (best split for mixed classes):

|N|masses only < 50×|< 27×|masses + deep voids < 50×|< 27×|< 10×|
|-|-|-|-|-|-|
|4|1|0|0|0|0|
|5|3|0|6|1|0|
|6|6|6|11|7|0|
|7|11|11|21|13|1|
|8|13|12|29|22|7|

- **Up to 4–5 objects the answer is close to unique.** 3 masses land in the same layout from most starts (121×). 4 masses land in a single basin (7 of 40 starts) at 35.6×.
- **From 6 objects up, almost every successful start is a different arrangement.** At N = 8, nearly every start below 27× is its own basin.
  - So beyond ~5 objects the details are weakly constrained (degenerate). Many unrelated layouts work about equally well.

## 3. What all good solutions share (`c_where_objects_go.png`, `features.json`)
Pooled over every distinct solution below 50× (1,666 masses and 784 deep voids):

- **Masses fall in four longitude bands:** ℓ ≈ 60–90°, 150–180°, 210–240° and 300–330°, about 80% of all masses.
  - These match the old four clumps' longitudes (75°, 162°, 228°, 315°), which form two mirror pairs about the 15°/195° axis.
  - This four-band backbone appears in every good solution, from the unique N = 4 fit to the many N = 8 variants.
- **Distance:** masses sit ~20 Gpc from us (10–90% range 16–28 Gpc); voids sit ~23 Gpc (17–31).
- **Out of the ring plane:** masses typically sit 15–25° out of the ring plane (median |b| 18°).
  - The in-plane spots at ℓ 162°/228° are inside the CoE interior, so the fit lifts those masses above or below it.
- **Voids sit between the mass bands:** near ℓ ≈ 0–30° (the side away from the CoE), 120–150° and 240–270°, and directly above and below the CoE direction (ℓ ≈ 195°, |b| ≈ 50°).
- **Migkas sector (ℓ 270–300°):** the fits mostly leave it empty. Only 3–4% of masses and 3–6% of voids land there. There is a void band at its edge, around ℓ ≈ 265–275°.
- **Size:** total cluster mass in the best fits is 3–7e52 kg, 1–2× the old four clumps' 3.1e52.
  - Size caps bind: 30% of masses sit at the 2e52 kg cap and 65% of voids at the 4 Gpc cap. Looser bounds would change the numbers, and earlier unbounded fits ran away to 1e54–1e58 kg.

## 4. Parameters vs constraints (when it becomes overfitting)
- **What has to be cancelled:** the residual is a pure cosine series about ℓ = 15.16° (all phases line up). The terms above the floor are m = 1–11 (m = 12 is 0.5×; everything past m = 12 adds 0.2×).
  - That makes ~11 numbers to cancel. Any layout that breaks the mirror symmetry also has to keep the ~11–12 sine terms near zero. In total: **12–24 effective constraints**.
- **What the fit can adjust:** each object has 4 numbers, so N = 3 → 12, N = 4 → 16, N = 6 → 24, N = 8 → 32.
- **Reading:**
  - N ≤ 4 is well constrained, and the unique 4-mass solution is informative.
  - N = 5–6 reaches the count of constraints.
  - N ≥ 6 has as many or more knobs than constraints. The gains from 6 → 8 objects (23× → 4.5×) come from fitting m = 7–12 detail, with many unrelated layouts giving the same quality. That is curve fitting, not information about real structure.

## Bottom line
- The robust result is a **backbone of 4 masses (~1–2e52 kg each, ~20 Gpc away, in mirror pairs at ℓ ≈ 75°/315° and 162°/228°, lifted out of the plane)**. It brings the leftover from 237× to ~36× and is essentially unique.
- Deep voids placed between those masses (away-from-CoE side, and above/below the CoE direction) push the leftover to 7–16×, but in many interchangeable ways.
- Shallow, self-compensated voids contribute nothing at any N.

## Files
- `ss_core.py`, `fit_run.py`, `analyze.py`, `features.py`, `report.py`
- `raw/*.json` (every start)
- `results.json` (per class: best, p10, median, fraction below thresholds, distinct basins)
- `features.json`
- Charts: `a_best_vs_N.png`, `b_distinct_below27.png`, `b_distinct_below50.png`, `c_where_objects_go.png`
- Logs: `log_fits.txt` (100-start batch), `log_fits2.txt` (40-start batch)

## 5. Keep adding until the floor (N > 8)
**Method:** same `ss_core` model, bounds and fitter. Each N starts from the best N−1 layout plus one new object (4 random spots for a new mass and 4 for a new void; masses-only track: 8 new-mass spots). All objects are then refit together. Stop rule: < 1×, N = 20, or 3 steps with < 2% gain each. Scripts: `extend.py`, `cold.py`, `polish.py`, `floor_check.py`, `report_ext.py`. Results: `extend_mixed.json`, `extend_masses.json`, `cold_*.json`.

**Masses + deep voids reach the floor at N = 11.**

|N|leftover|added|largest remaining modes|at cap (mass / void)|
|-|-|-|-|-|
|8|4.46×|start|m = 10, 9, 12|1 / 3|
|9|2.96×|mass|m = 10, 9, 12|1 / 2|
|10|1.26×|mass (void 1.30×)|m = 10, 12, 11|0 / 2|
|11|**0.58×**|void (mass 0.59×)|m = 12, 11, 10 (all ≤ 0.09×)|0 / 3|

**Masses only never get there.** 13× (N = 8) → 11.2, 9.7, 8.4, 7.7, 7.1, 6.3, 5.7, 5.3, 4.9, 4.6, 4.3 → **4.1× at N = 20**. Gains are ~5–10% per added mass and still falling, so reaching 1× would take far more than 20. The stuck modes are m = 10, 7 and 11/12, and no mass sits at the cap. Without voids, single masses can't make the m ≈ 7 and 10 dips.

**The floor solution:** 6 masses + 5 deep voids, 0.58×. Every m = 1–12 term is ≤ 0.09×.
- **Mirror symmetry:** it is almost exactly mirror-symmetric about the ℓ 15°/195° axis. There are 5 mirror pairs plus one void on the axis.
  - Mass pairs: ℓ 160°/230° (1.84e52 kg, 21.4 Gpc), 175°/215° (6.4e51, 16.9 Gpc) and 85°/305° (9.1e51, 17.7 Gpc, |b| 35°).
  - Void pairs: ℓ 351°/39° (Rv 4.0, 32 Gpc, in plane) and 123°/267° (Rv 3.5, 17.8 Gpc, |b| 32°). The axis void is at ℓ 15° (Rv 4.0, 30.8 Gpc).
- **Same backbone:** the four-band mass pattern seen at N ≤ 8 is still there.
  - The pairs at 160°/230° and 175°/215° sit right on the CoE-interior edge (44.5–44.7 Gpc from the CoE, against the 43.9 Gpc limit plus core size). So the exclusion boundary is binding.
- **Mass balance:** clusters total **6.77e52 kg** and the void deficit is **−6.76e52 kg**. The net is +1.4e50 kg (0.2%), so the layout balances its own mass.
- **Caps:** 3 of 5 voids are at the 4 Gpc cap. No mass is at the 2e52 kg cap.

**Parameters vs constraints at the floor:**
- 11 objects × 4 = **44 free numbers**.
- The residual has ~11 cosine terms above the floor (m = 1–11), or ~24 numbers counting the sine terms that have to stay at zero.
- The symmetric layout really uses ~22 independent numbers (5 pairs × 4 + the axis void's distance and size).
- So the fit has roughly 2× (by the symmetric count) to 4× (by raw count) more freedom than the ~11 numbers it must cancel. **Reaching 1× at N = 11 is a fit, not a prediction.** The informative part is the backbone found at N ≤ 5.

**Degeneracy (cold random starts):**
- 6 masses + 5 voids: 16 cold starts, extra-polished. Best 1.07×, and **0 of 16 below 1×**. Nine are under 3×, all in different layouts.
- 6 masses + 6 voids: 16 cold starts. **1 of 16 below 1×** (0.72×).
- **Below 1×:** 2 distinct solutions in total (the warm-path 0.58× and this 0.72×).
- **Below 2×:** 15 solutions in 13 distinct layouts.
- So near the floor there are many separate, nearly equal layouts. The cold starts are hard to push under 1× from scratch; the warm path, which keeps the backbone, gets there more easily.

**Cuts:** cold starts were 16 per N, and only at N = 11 and 12. The new object got 4 random spots per type per step, not more.

Charts: `a_best_vs_N_logy_extended.png`, `d_floor_solution_map.png`
