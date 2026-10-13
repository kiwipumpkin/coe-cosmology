# Rerun with the CoE core fully open

Chris: "We specifically asked to eliminate that restriction. Please run them once again without the core excluded."

## What changed
- **Nothing is excluded at the core.** There is no wall and no smooth-density mass cap. Masses are 1e49–2e52 kg, deep voids (δ −0.7) have radius 0.5–4 Gpc, and giant partners are 1–15 Gpc. These are the same bounds used everywhere else.
- **Reach:** along every direction a body can go out to max(43 Gpc, the far edge of the core). Toward the CoE (ℓ 195°) the core runs from 14.3 to **102 Gpc** from us, so the whole core, far side included, is reachable. For partners the old 60 Gpc reach is kept away from the core.
- **Kernels, compensation and other rules are unchanged:** `ss_core` (outside, no-mirror), `ao_core` (all open, no density limit), `ep_core` (partners, strict net mass zero, outer shell edge kept for partners). Code: `co_core.py`.
- **Starts (214 fits):** warm (from the earlier best fits), cold (uniform over the whole open reach), and **core-seeded**, where 1/3, 1/2 or all of the bodies start at random points inside the core.
- **Extra check:** fits with *every* body forced inside the core (`run_coreonly.py`, 6 starts per mix).

## Best leftover (× floor) — `best_vs_N_core_open_logy.png`
| N | outside: walled → core open | no mirror: walled → open | giant partners: walled → open | all open: earlier → open | all bodies in core |
|-|-|-|-|-|-|
| 4 | 35.6 → 33.7 | | | 7.27 → 7.27 | 172 |
| 6 | 15.9 → 15.5 | | 10.4 → 8.7 | 1.41 → 1.40 | 172 |
| 8 | 4.46 → 2.97 | | | **0.18 → 0.18** | 154 |
| 10 | 1.26 → **0.46** | 2.30 → 1.24 | 1.23 → **0.17** | | 151 |
| 11 | 0.58 → 0.24 | 1.25 → **0.37** | — → 0.080 | | |
| 12 | 0.72 → 0.069 | 0.73 → 0.34 | — → 0.49* | | 143 |

Also: 3 partners alone, 106 → 40. *Undersampled: the N = 12 partner fit, the outside 5+4 mix at N = 9 (14.7×, 2.96× walled with 6+3), and all-open 7 masses (13.7×, never run before) had no good warm starts. They are in `report_out.txt`, but N = 9 and all-open N = 7 are left off the chart.

- **First N at or below 1×:** outside 10 (was 11), no-mirror 11 (was 12), partners 10 (never reached before), all-open 8 (unchanged).
- **Almost all of this gain was already there with the wall off at the old 43 Gpc reach** (`../solution_space/no_core_wall`: 0.46 / 0.25 / 0.074 at N = 10 / 11 / 12). Opening the deep core added nothing: 0.455 / 0.241 / 0.069.

## Where the fitted bodies go
Pooled over all fits below 2×:
- **Outside:** 26% of bodies are in the core (24 of 78 masses, 14 of 66 voids). Their median depth is **3 Gpc** and their median distance from us is 24 Gpc, at ℓ ≈ 196°. Only 4 are deeper than 10 Gpc, and only 3 are farther than 50 Gpc from us.
- **No mirror:** 33% in the core, median depth 4 Gpc, all within 32 Gpc of us.
- **Partners:** 28% in the core, median depth 6 Gpc, mostly within 30 Gpc. In the best fit, partner P2 is an overdense sphere (R 11 Gpc, +3.1e52 kg) sitting 9.5 Gpc inside the core edge at ℓ 225°.
- **All open:** 4% in the core. The best fits put all 8 bodies 5–7 Gpc from us, inside our sphere.
- **Core-seeded starts:** they move their bodies back out to the near edge, or leave them parked uselessly. They win only 4 of 24 classes and score a median 4–24× worse than warm starts.
- **Core-only fits:** they never get below **143×** (12 bodies), against 237× at the start. The fitted bodies hug the near edge of the core (0–4 Gpc deep, 16–37 Gpc from us). Extra masses get parked on the far edge at 102 Gpc, where they do nothing.

Best layouts (side views, with the leftover around the ring on its own scale): `side_view_out.png` (12 bodies, 0.069×, 2 in core), `side_view_asym.png` (12, 0.34×, 3), `side_view_ep.png` (11, 0.080×, 3), `side_view_ao.png` (8, 0.18×, 0). Body-by-body tables are in `report_out.txt`.

## Do core bodies help? Only at the near edge, and the deep core doesn't
- **The near edge matters.** It sits 15–30 Gpc from us toward ℓ 160–230°, which is exactly where the old 162°/228° mass bands and the ℓ 195° axis void want to be. That's why lifting the wall helped: 0.72× → 0.07× at 12 bodies. In the best outside fit, the 2 edge bodies are needed. Removing them without refitting sends the leftover back to 197×.
- **The deep core is useless, for a simple reason.** A body's pull on the ring depends on how much the pull changes across the ring. A body far away pulls almost equally on every part of the ring, so it leaves almost no pattern, and what it does leave is just a gentle two-lobed (m = 2) stretch. The m = 3–5 wiggles we need to cancel (200, 185, 173×) die off much faster with distance. One 2e52 kg mass on the CoE line gives:

| distance from us | m = 2 | m = 3 | m = 4 | m = 5 |
|-|-|-|-|-|
| 20 Gpc (edge) | 1036 | 862 | 691 | 541 |
| 30 Gpc | 205 | 119 | 65 | 35 |
| 45 Gpc | 53 | 21 | 8 | 3 |
| 60 Gpc (CoE) | 22 | 6 | 1.8 | 0.5 |
| 100 Gpc (far edge) | 4.5 | 0.8 | 0.1 | 0.0 |

So a body at the CoE is about 100–1000× weaker in the modes that matter than the same body at the edge. Matching the target from there would take masses far above the 2e52 cap, and those would bring a huge m = 2 that would then need cancelling too. The optimizer was given the whole core and seeded inside it, and it still chose the edge.

## Caveats
- Starts were limited by time (8–12 per mix). Some classes (marked * above) are clearly undersampled.
- The far-side reach (core exit at 102 Gpc) is not limited by the outer shell edge for ordinary bodies, the same as the old 43 Gpc reach.
- Same Newtonian kernel as before: no ISW/Rees–Sciama or lensing.

Files: `co_core.py`, `run_co.py`, `run_coreonly.py`, `report_co.py`, `raw/`, `raw_coreonly/`, `coreonly_best.json`, `best.json`, `report_out.txt`, `log_run.txt`, `log_coreonly.txt`. The brief and notes/ are untouched.
