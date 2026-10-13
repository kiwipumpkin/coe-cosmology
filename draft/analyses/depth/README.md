# Tidal influence depth ("multipole horizon")

Question (Chris): beyond some shell of clusters and voids around us, does the distribution stop mattering? How deep is that shell?

## Code
- `fits_cutoff.py` (A): refits k placed objects with distances limited to [15.0 Gpc, R_cut]. Usage: `python3 fits_cutoff.py <k> <signed> <coe> <tag>` → `fits_<tag>.json`, `log_<tag>.txt`.
  - It imports `../shell_ensemble/void_mix/followup/fit_masses_voids.py` unchanged: same model, net-mass compensation, size cap (void at most empty, δ=−1) and no-overlap rule.
  - Only the distance bound changes. `coe=1` also keeps every object at least D−d = 43.9 Gpc from the CoE, the same region as the random web. That is what allows going past 42.9 Gpc.
  - Warm starts: the previous R_cut solution plus the unrestricted fit, plus 120 random starts.
- `web_shells.py` (B, C) → `web_shells.json`, `cumulative_curves.json`.
  - B1 regenerates the void_mix k=0, N=200, medium, 50%-void layouts exactly (max difference vs the stored leftover is 3e-17), then splits each layout into distance shells.
  - B2 uses the same object population (uniform in volume, CoE interior excluded) extended to 400 Gpc, to get the share of the effect inside R for a web that doesn't stop at 42.9 Gpc.
  - C covers the analytic falloff and the measured slopes.
- `report.py` → `results.json`, `summary_table.md`, `depth_cutoff.png`, `effect_vs_distance.png`.

## Definitions
- "Effect" of the random web = the power it adds to the leftover. That is the mean square of the web's pattern; the earlier check gave E[leftover²] = start² + E[web²].
  - Because positions and signs are random, this adds object by object, so the "share from inside R" is well defined.
  - Mode groups use the same power, restricted to m=1–2 or m=3–5.
- Depth R90 / R95 = the distance from us containing 90% / 95% of that effect.
- For placed objects: R_cut where the best leftover is within 5% of its final value, and where 90% / 95% of the total improvement (from 237×) is reached.

## What sets the depth
- **The field model:** the tidal pull across our visible sphere, with the uniform pull (dipole, l=1) removed. A mass at distance r produces ring mode m with amplitude ∝ (R_vis/r)^(m+1) for m≥2 and ∝ (R_vis/r)^4 for m=1. The m=1 mode only appears at l=3 because l=1 is removed.
  - Measured slopes at 60–180 Gpc are −4.01, −2.95, −3.93, −4.91 and −5.88 for m=1–5, which match. Close in (15–30 Gpc) the falloff is about one power steeper.
- **The population:** a web that is uniform in volume grows as r², so the effect per distance falls as r^(−2m) for m≥2. Starting at r0, the share inside R is 1 − (r0/R)^(2m−1).
  - m=2 is the slowest and gives R90 = 2.15 r0. The far-field formula is an upper bound, because the real near field is steeper.
- **The gap:** the decisive number is how close structure gets to our 14.3 Gpc sphere. Near the sightline ring a single object's effect grows like 1/gap², so the effect piles up right next to the sphere.
  - With object centres allowed right at 14.3 Gpc, 99.7% of the random web's effect comes from 14.3–20 Gpc.
  - Requiring objects to sit wholly outside the sphere (start at about 15.0 Gpc for medium-mass objects) gives the headline numbers.
- **Placed objects:** these can trade distance for mass (pull ∝ M/r³), so their depth is softer. For 4 objects the fit flattens at 26 Gpc. For 3, the best fit uses one very deep, far void (1.4e53 kg at 50 Gpc) for the m=2 part, and only flattens near 42 Gpc.
  - If the CoE interior is allowed, the 4-object "masses or voids" fit reaches 52× using a void inside D−d (the earlier followup result). With the interior excluded, the best 4-object fit is all clusters, at 73×.
