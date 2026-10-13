# Rules for random void + cluster networks (null ensemble)

**What this is.** Rules for filling a region of space with random voids and galaxy clusters that look like the real
cosmic web. They're meant to be the "nothing special" comparison set (null ensemble) for the CMB statistics.
The links only join voids to clusters (bipartite). A void links to every cluster on its boundary, and a cluster
links to every void it touches. There are no void–void or cluster–cluster links.

**How the generator builds the web (picture from Icke & van de Weygaert 1987 and van de Weygaert & Icke 1989).**
Voids are the cells of a foam. Walls are where two cells meet, filaments are where three meet, and clusters sit at
the corners (nodes) where four or more meet. The generator:
1. drops void seeds,
2. turns each seed's cell into a void,
3. puts web nodes at the cell corners,
4. makes the most prominent nodes into clusters, with masses taken from the measured cluster mass function, and
5. links each cluster to the voids around it.

"ours" means we found no published measurement, so the value is our own plausible choice and can be varied.
"derived" means the number follows from a cited result by a short calculation, which is shown.

## Parameter table

| # | parameter | value / range used | source |
|---|---|---|---|
| a | **clusters bordering each void**: distribution | mean **8** (range we'd allow: 4–27). The spread comes out of the foam: roughly bell-shaped, about 1–18 (see chart). | **ours**. We found no published count. Upper limit: a random (Poisson–Voronoi) foam cell has on average 27.07 corners (Meijering 1953, *Philips Res. Rep.* 8, 270; pre-arXiv), so 27 would mean every corner holds a cluster. The mean has to obey the counting rule (clusters per void) × (number of voids) = (voids per cluster) × (number of clusters). |
| b | **voids touching each cluster**: distribution | mean **4.5** (allowed 4–5). Distribution: 4 (66%), 5 (23%), 6 (8%), 7+ (3%). Rises with cluster mass. | **derived** from Codis, Pogosyan & Pichon 2018, **arXiv:1803.11477**: a peak has about 4 filaments locally and about 6.1 on average globally in 3D, and this rises with peak height (about 3 up to more than 10). Around a node where every filament borders 3 walls, Euler's rule on a small sphere gives voids = 2 + filaments/2. That makes 4 filaments → 4 voids and 6 → 5. The rise with mass agrees with Aragón-Calvo, van de Weygaert & Jones 2010, **arXiv:1007.0742** (about 2 filaments at 10^14 Msun/h, about 5 at ≥10^15 Msun/h), and with observations by Sarron et al. 2019, **arXiv:1903.02879** and Darragh-Ford et al. 2019, **arXiv:1904.11859**. |
| c | **cluster distance from void centre**, in void radii (Rv = equal-volume radius, as in VIDE) | not set directly; the foam fixes it. Measured: median **1.31 Rv** (68% range 1.15–1.49). | Lower bound: the overdense ridge around a void peaks at **1.09 Rv**. That's our calculation from the Hamaus, Sutter & Wandelt 2014 profile (**arXiv:1403.5499**) at the compensated values rs/Rv = 0.91, α = 2, β = 9.5. Their walls peak "slightly outside Rv", and clusters sit at the corners, which are farther out than the wall midpoints. Rv definition: Sutter et al. 2015 VIDE, **arXiv:1406.1191**; ZOBOV, Neyrinck 2008, **arXiv:0712.3049**. |
| d | **angular spread of a void's clusters** on its boundary | the foam fixes it. Measured: a void's clusters are spread about as evenly as random directions would be (median lopsidedness 0.34, against 0.32 for random directions with the same counts). | **ours** (no measurement found) |
| e | **cluster mass vs size of the bordering voids** | rank nodes by the size of their empty sphere (bigger empty corner means bigger surrounding voids) with ±0.5 (natural-log) random scatter, then give masses in rank order. Measured slope: M ∝ (mean bordering Rv)^**0.25**. | **ours** (no measurement found). The rising trend matches row (b): bigger clusters touch more voids. |
| f | **compensation** (void's missing mass vs mass in its ridge) | default "lcdm": each void cell plus its ridge totals the cosmic mean up to a random ±σ_LCDM(cell radius) (about ±1% for 480 Mpc cells). Variants: "exact" (all cells exactly compensated); "hamaus" (bigger voids under-compensated, f = 1 − 0.05 ln(Rv/25.1 Mpc)). | Literature: voids of about 17.6 Mpc/h are exactly compensated, bigger ones under- and smaller ones over-compensated (Hamaus et al. 2014, **arXiv:1403.5499**; void-in-void / void-in-cloud: Sheth & van de Weygaert 2004, **astro-ph/0311260**). Lensing measures about −0.4 inside Rv, partly compensated, with a slow approach to the mean (Clampitt & Jain 2015, **arXiv:1404.1834**). Choosing σ_LCDM for the cell-to-cell scatter, and the 0.05 slope, is **ours**. |
| f2 | mass budget of the web | voids fill **77%** of the volume at mean density 1+δ = **0.2**. The remaining mass (walls, filaments, smaller haloes) sits on each void's ridge. | Cautun et al. 2014, **arXiv:1401.7866** (NEXUS+ table: clusters 11% of mass and <0.1% of volume; filaments 50% / 6%; walls 24% / 18%; voids 15% / 77%; void mean 1+δ = 0.2) |
| g1 | **cluster mass function** | Tinker et al. 2008, Δ = 200 × mean density, z = 0 (A = 0.186, a = 1.47, b = 2.57, c = 1.19), on the same linear σ(R) as density_vs_scale. Masses above Mmin (demo: 3×10^15 and 2×10^15 Msun) | Tinker et al. 2008, **arXiv:0803.2706** |
| g2 | **void size function** | the foam sets it: Rv ≈ (3 V_cell / 4π)^(1/3), median **479 Mpc** at Mmin = 3×10^15 (cells around the rarest clusters, so these are supervoid-scale voids). Number of voids = (clusters × 4.5) / 8. | **ours** at these scales. Dark-matter void abundance models (Sheth & van de Weygaert 2004, **astro-ph/0311260**; Jennings, Li & Hu 2013 "Vdn", **arXiv:1304.6087**, coded in lss.py) cover voids of about 1–15 Mpc/h in the dark matter. Jennings et al. note that voids traced by haloes are larger and more numerous at r > 10 Mpc/h, so those models don't apply directly to voids bounded by >10^15 Msun clusters. |
| h | non-overlap | void cells never overlap (they tile space). The drawn void "core" spheres are shrunk so neighbours never overlap: each pair splits its separation in proportion to their sizes. Seeds have a soft minimum spacing of 0.7 × mean spacing, which makes the foam more regular than pure random. | tiling: Voronoi construction. Minimum spacing 0.7: **ours** |
| i | node merging | foam corners closer than 0.15 × mean seed spacing are merged into one node, which gives nodes with 5, 6 or more voids. Set so that the mean of (b) is 4.5. | **ours** (calibrated to row b) |
| j | region and exclusions | ball of radius 20 Gpc (our 14.3 Gpc sphere plus the shell out to 18 Gpc and a margin). Nothing within 0.5 Gpc of us, the same gap as the joint fits. No lightcone evolution: everything is at z = 0. | **ours** |
| k | **links never cross** (Chris's rule) | Links are tubes of radius **24 Mpc** (0.05 × median void radius). Two links that don't share an end must stay more than 2 tube radii apart. No link may pass through another node's body: a void's core sphere, or a cluster's R200m sphere. That covers cluster–cluster links crossing voids. Candidates are tried shortest first. If one crosses a link already kept, it is dropped, so of two crossing links the longer one goes. One that passes through a body is dropped outright. | **ours** (rule from Chris; tube radius ours, tested at 10 and 48 Mpc). Colberg, Krughoff & Connolly 2005 (**astro-ph/0406665**) note that requiring cluster–cluster connections not to cross voids would raise the fraction that are real filaments. |

## What the check shows (`stats_check.png`, `stats_check_3e+15.json`; one network, Mmin = 3×10^15 Msun)
- (a) clusters per void: mean 7.89 (target 8)
- (b) voids per cluster: mean 4.50 (target 4–5); 4/5/6/7+ = 66/23/8/3%; mean rises with mass from 4.5 to 5.2
- (c) cluster–void distance: median 1.31 Rv, always outside the 1.09 Rv ridge peak
- (d) a void's clusters are spread about as evenly as random directions with the same counts
- (e) slope 0.25 (free; no target)
- (f) net contrast per cell: rms 0.0103, against LCDM σ 0.0099 at the cell size
- (f2) void volume fraction 0.770 (target 0.77)
- (g1) the generated cluster mass function matches Tinker08 at every mass; the mass fraction in clusters is 3.36e-4 against 3.37e-4 from Tinker

## Files
- `lss.py`: void abundance (SvdW + Vdn) and cluster mass function (Tinker08) on the shared σ(R)
- `vmn_gen.py`: the generator (`generate(**params)`, `kernel_bodies(g)` gives point masses for the CMB kernel)
- `check_stats.py` produces `stats_check.png`; `side_view.py` produces `network_side_view.png`
- `planar.py`: the no-crossing filter (`planarize(g, cc, cc2)`). `census.py` runs it and writes `census_3e+15.json` with before/after numbers. `census_plots.py` makes `census_counts.png`, `census_lengths.png` and `census_angles.png`. `planar_view.py` makes `side_view_noncrossing.png`
- `demo_ensemble.py` / `demo_plot.py` produce `results/demo_*.json`, `demo_leftover_hist.png` and `demo_summary.json`

## Caveats (plain)
- Clusters are the rarest ones only (above 2–3×10^15 Msun), so the voids are the big cells between them. Going down to
  ordinary 10^14 Msun clusters would need about 10^8 objects. The generator handles any Mmin, but below about
  1×10^15 Msun the whole 20 Gpc ball needs splitting into tiles.
- Cells are placed independently (a foam with no extra large-scale clustering). Large-scale LCDM power beyond the cell
  size is therefore missing. The demo spread is a lower bound on how much random structure can move the CMB.
- Masses are point masses at cell centroids and nodes. A ridge spread over the cell boundary acts like a point at the
  centroid only when seen from far away. Cells next to the 0.5 Gpc gap are the main exception.

## Demo null ensemble through the CMB kernel (`demo_leftover_hist.png`, `demo_summary.json`)
Leftover = rms(CMB residual − network pattern) / floor. With no structure at all it is 237×; the target is 1×.

| ensemble | n | leftover median (min–max) | rms of the network's own pattern (median) |
|---|---|---|---|
| clusters >3e15, cells compensated to LCDM σ (default) | 50 | 245× (235–421) | 57× (inside our sphere 38×, outside 31×) |
| clusters >2e15, same | 6 | 249× (244–256) | 73× |
| clusters >3e15, every cell exactly compensated | 10 | 237.3× (237.1–237.5) | 2.7× |
| clusters >3e15, Hamaus trend extrapolated | 10 | 478× (420–741) | 427× |

Reading: random networks that follow the rules have no tendency to cancel the residual. They add a pattern that has
nothing to do with it, so the leftover stays at 237× or gets worse. 4 of the 50 default networks came out a little
below 237× (best 235×). How big the added pattern is depends almost entirely on how well each void is compensated,
not on the clusters. Exact compensation gives about 3×, LCDM-sized scatter about 60×, and stretching the Hamaus
under-compensation out to 500 Mpc voids gives about 400×. None of the random networks gets anywhere near 1×.

## Link census with the no-crossing rule (`census_3e+15.json`, one network, Mmin = 3×10^15 Msun)
We searched for crossings inside 7 Gpc of us. We counted links whose midpoint is inside 4.5 Gpc: 1522 clusters and 819 voids. The table shows how many links the rule removes, by type, at a tube radius of 24 Mpc:

| link type | candidates | removed: crossing | removed: through a body | removed % | removed % at 10 / 48 Mpc tubes |
|---|---|---|---|---|---|
| cluster–cluster, direct (foam edge) | 1677 | 12 | 1 | 0.8% | 0.5% / 2.1% |
| cluster–cluster, via one empty node | 2690 | 220 | 159 (through a void core) | 14.1% | 9.7% / 22.9% |
| cluster–void | 6879 | 21 | 0 | 0.3% | 0.1% / 1.7% |
| void–void (shared wall) | 6187 | 1167 | 15 | 19.1% | 9.0% / 30.4% |

What changes, before → after:
- filaments per cluster (both kinds): 5.71 → **5.20**. Direct only: 2.19 → 2.17
- voids per cluster: 4.51 → 4.49. Clusters per void: 8.47 → 8.44. The bipartite links are almost untouched.
- neighbour voids per void: 15.22 → **12.27**. This is the big change: a void–void link usually passes close to one of the cluster–void links that ring their shared wall.
- median lengths: cc direct 378 → 375 Mpc; cc via one node 762 → 717 Mpc; cluster–void 646 → 645 Mpc; void–void 980 → 923 Mpc (in void radii: 2.03 → 1.92). The rule removes the long links first.
- angles at nodes barely move (filaments at a cluster: median 98° → 101°).

The CMB kernel only sees the point masses of voids and clusters, not the links. So this rule does not change the demo leftover histogram. It only changes the network statistics above.
