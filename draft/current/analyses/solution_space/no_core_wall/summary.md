# Removing the hard wall at the CoE interior

Chris: "The density fx is smooth." The outside fits had treated the CoE interior (within 43.9 Gpc of the CoE) as a forbidden zone. This rerun removes that wall.

## Density used
The model's own exp-lock density (`local_web/web_core.rho_gpc`):
- **ρ(r) = 2.9×10⁻²⁷ kg/m³ × exp(−(r − 58.2 Gpc) / 922.65 Gpc)**, where r is the distance from the CoE.
- It is smooth everywhere and has no edge at 43.9 Gpc. Across the fit region it varies by only a few percent: ρ(43.9) = 1.016 ρ̄ and ρ(0) = 1.065 ρ̄.

## Variants (everything else unchanged)
Unchanged: the 14.3 Gpc sphere, the 43 Gpc reach, the 1e49–2e52 kg mass range, 0.5–4 Gpc deep voids at δ −0.7, and the smooth net-mass compensation.
- **Old:** the hard wall (`../` results).
- **Smooth-density limit:** the wall is removed. Each mass is limited by what the local smooth density can supply:
  - **Cap:** M ≤ 2e52 kg × ρ(x)/ρ̄.
  - **Supply sphere:** a mass draws its material from a sphere of volume M/ρ(x), about 3.8 Gpc radius for 2e52 kg.
  - **No double counting:** supply spheres and void spheres may not overlap each other or reach into our visible sphere.
  - **Voids:** their deficit is 0.7·ρ(x)·V, which already follows the local density.
- **Fully open:** the wall is removed and nothing replaces it, beyond the usual caps and cores not overlapping.
- **Starts:** 12–32 per N (warm starts from the old fits plus cold starts). The no-mirror rule from `../asym` was also applied for N = 11 and 12 under the smooth limit, with 10 starts each.

## Results (`a_best_vs_N_no_wall_logy.png`)
|N|old (wall)|smooth limit|fully open|smooth + no mirror|
|-|-|-|-|-|
|4|35.6|33.9|33.7| |
|6|15.9|14.2|14.2| |
|8|4.46|2.91|2.97| |
|10|1.26|**0.73**|0.46| |
|11|0.58|0.17|0.25|0.54|
|12|0.72|0.094|0.074|0.34 (no-mirror old: 0.73)|

- **First below 1×:** at N = 10 (6 masses + 4 voids), one object fewer than with the wall.
- **Many solutions below 1×:** 14 of 16 starts at N = 11 and 15 of 16 at N = 12.
- **N ≤ 6 barely changes:** the four-mass layout (33.9×) is the same one as before. Its two masses at ℓ 162° and 228° now sit just 0.4 Gpc inside the old wall line, so the wall was not what limited it.
- **Smooth limit vs fully open:** almost the same. The density varies only a few percent, and the supply-sphere rule rarely binds.

## Where things go without the wall
Pooled over every fit below 2×:
- **Masses:** about 31% now sit inside the old wall line. They go only a little way in (median 2.4 Gpc, at most 8 Gpc), at ℓ ≈ 160–230° about 20–25 Gpc from us.
- **Voids:** 17% sit inside, with a median depth of 6 Gpc and up to 19.5 Gpc. A typical example is a 4 Gpc void on the us–CoE line (ℓ 195°, b 0°) about 12 Gpc inside the old wall.
- **The masses don't pile up deep in the core.** The fit moves a few masses slightly past the old line and puts a big void right on the CoE axis.

## Map: smooth limit, first N below 1× (N = 10, 0.73×; `d_no_wall_smooth_map.png`)
- **Masses:**
  - ℓ 55° and 83°: 33–35 Gpc out, both at the cap.
  - ℓ 164° and 224°: inside the old wall by 4.4 and 2.5 Gpc.
  - ℓ 258°: a small one.
  - ℓ 323°, b +39°.
- **Voids:**
  - ℓ 195° on the axis, 11.9 Gpc inside the old wall.
  - ℓ 128°, 258° and 280°.
- **Mass balance:** 6.13e52 kg in masses against −6.11e52 kg in voids, which agree to 0.3%.
- **Objects at limits:** 2 masses at the density-scaled cap and 3 of 4 voids at the 4 Gpc cap. Nothing is at the 43 Gpc reach.

Files:
- Code: `ncw_core.py`, `run_ncw.py`, `report_ncw.py`
- Results: `raw_smooth/`, `raw_open/`, `raw_smooth_asym/` (`raw_smooth_asym_part1/` holds 4 duplicate starts from an interrupted launch), `ncw_best.json`, `report_out.txt`
