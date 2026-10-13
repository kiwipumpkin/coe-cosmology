# Everything open: objects anywhere (inside our sphere, in the shell, in the old core)

Chris: the no-wall run should also drop the restriction on our visible sphere.

## What's allowed
- **Region:** objects 0.5–43 Gpc from us in any direction. There is no CoE wall and no visible-sphere exclusion.
- **Near us:** no core or void sphere may come within 0.5 Gpc of us.
- **The CMB ring (14.3 Gpc): straddling is forbidden, not split.** Each core or void sphere must lie wholly inside (d + R ≤ 14.2) or wholly outside (d − R ≥ 14.4).
- **Kernel:** `ss_core.unit_field` (= `field3`) is used for every object. It is plain Newtonian (pull at the ring point minus pull on us, projected on the line of sight), and `inside_visible` verified it for interior points. So the inside and outside kernels are the same function.
- **Density:** ρ(r) = 2.9e-27 kg/m³ × exp(−(r − 58.2)/922.65), with r in Gpc from the CoE.
  - **Smooth-density limit:** a mass is capped at M ≤ 2e52 kg × ρ(x)/ρ̄. Each mass draws from a supply sphere of volume M/ρ(x), and supply spheres and void spheres may not overlap. Supply spheres may cross the ring and enclose us; only the physical cores and voids are restricted.
  - **No density limit:** only the usual caps, and cores may not overlap.
- **Compensation (`inside_visible` convention, used for every fit here):** the net mass is removed smoothly over our visible sphere, following the exp-lock density shape (`inside_core.FLOC`, the analytic gradient dipole). A uniform sphere contributes exactly zero.
  - The comparison curves keep their own conventions: inside-only also uses this local one, while outside-only (no wall) uses the old smooth compensation from the surround.
- **Starts:** 12 per mix (8 for the no-density-limit run). Starts were warm from the inside-only best, warm from the outside-no-wall best, mixed (each object taken at random from one or the other), or cold random.
- **Mixes:** masses only, plus mixed splits at N = 4, 6, 8, 10, 11, 12.

## Result (`a_best_vs_N_all_open_logy.png`)
|N|inside only|outside only (no wall)|all open, smooth limit|all open, no density limit|
|-|-|-|-|-|
|4|7.27|33.9|7.28|7.27|
|6|2.20|14.2|1.41|1.41|
|8|0.18|2.91|**0.29**|0.18|
|10|0.060|0.73|0.15|0.059|
|11|0.051|0.17|0.090|0.048|
|12|0.025|0.094|0.073|0.025|

- **First N below 1×:** N = 8 (4 masses + 4 voids, 0.29×), the same as inside-only. N = 7 was not run.
- **No-mirror rule at N = 8:** the best is **3.66×** (4 + 4) and 6.1× (8 masses), from 12 starts each, so it no longer reaches 1× at N = 8.
- **With no density limit,** all-open reproduces inside-only exactly. The best fits put every object inside our sphere.
- **With the smooth limit,** all-open does worse than inside-only from N = 8 up. Supply spheres of about 3–4 Gpc can't all fit near us without overlapping, so the packing is limited.
  - At N = 11 and 12 the best smooth fits switch to the all-outside layout (0.090×, 0.073×). That beats the earlier outside-only fits run with the old compensation (0.17×, 0.094×).

## Where objects go when everything is open
- **Best fits are all one side:** every best fit puts all its objects inside our sphere (N ≤ 10) or all of them outside (N = 11, 12, smooth limit). None of the best fits is a mixture.
- **Pooled over all smooth-limit fits below 2×:** 61% of masses and 42% of voids are inside our sphere.
  - **Inside:** median distance from us is 6.9 Gpc for masses and 8.1 Gpc for voids.
  - **Outside:** median distance is 27 Gpc.
  - **Old core:** only 11% of objects sit inside the old wall line.
- **Starts that mix both sides do worse:** mixed and cold starts end with a mixture (about 30–90% of objects inside), and score worse than the pure inside basin (e.g. N = 8: 0.65× mixed, 1.57× cold, 0.29× inside).

## Map: N = 8, smooth limit, 0.29× (`d_all_open_map.png`)
- **Positions:** all 8 objects are inside our sphere, 5–10 Gpc from us, spread around the sky.
  - Masses: ℓ 16°/+37°, 128°/−35°, 195°/+30°, 262°/+23°.
  - Voids: ℓ 80°, 165°, 226°, 308°.
- **Mass balance:** masses total 6.78e52 kg and voids −3.36e52 kg. The net +3.4e52 kg (3% of the visible sphere's 1.04e54 kg) is removed by the local compensation.
- **Objects at limits:** 1 mass at the density-scaled cap. No voids at the 4 Gpc cap, nothing at the ring margin, the near-us margin or the 43 Gpc reach.
- **What's left:** only m = 7–12, all at or below 0.22×.

Caveat (from `inside_visible`): for objects between us and the CMB, the model leaves out the ISW/Rees–Sciama effect, lensing and potential-difference effects, which are probably not small at 1e52 kg and 5–10 Gpc.

Files:
- Code: `ao_core.py`, `run_ao.py`, `report_ao.py`
- Results: `raw_smooth/`, `raw_open/`, `raw_smooth_asym/`, `report_out.txt`
