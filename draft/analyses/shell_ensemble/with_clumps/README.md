# Random web on top of k = 0–4 placed clumps (Chris: "place one through four of the required masses, mindful of the density gradient")

Scripts: `baselines.py` (clump patterns, best subsets, compensation, refit checks → `baselines.json`, `baseline_patterns.npz`),
`ensemble_with_clumps.py` (run: `python3 ensemble_with_clumps.py 1000` → `counts.json`, `raw_*.npz`), `report.py` (→ `summary_table.md`, figures).
Random web: identical to `../surround` grad_bare (same region, exp-lock gradient, N, mass levels, zero net contrast, seeds). k=0 reproduces that run exactly.

- **Clumps:** the fitted four from `scan.json` sym N=4 (83.4× fit). Masses are 7.52e51, 7.52e51, 7.79e51 and 7.79e51 kg. All four sit 24.5 Gpc from us, in the ring plane.
  - Clumps 0 and 1 (ℓ 162°/228°) are 39.9 Gpc from the CoE, which is inside D−d.
  - Clumps 2 and 3 (ℓ 75°/315°) are 73.6 Gpc from the CoE, just beyond D+d. None lies inside the 43.9–72.5 Gpc shell itself.
- **Subsets:** the best subset of the fitted four, keeping their fitted masses and positions. The refits listed below are for information only.
- **Compensation (adopted):** each placed clump's mass is removed from the surround region (14.3–42.9 Gpc from us, CoE distance ≥ D−d) in proportion to ρ0 e^(−r/L).
  This scales the smooth profile by (1−ε), with ε ≈ 0.1%, so its shape and gradient are kept and the net added mass is zero.
  The deficit's tidal pattern was computed with the same field3, using 400k sample points.
- **Bracket, local compensation:** a deficit concentric with each clump (radius 2.8–9 Gpc, not reaching our 14.3 Gpc sphere) cancels that clump's pattern. Shell theorem check: 170.7× becomes 0.1–0.7×, which is MC noise.
  So fully local compensation returns every baseline to 237×.
- **Overlaps:** random objects whose sphere would touch a placed clump's δ=+1 sphere (2.8 Gpc) are redrawn. This affects about 0.1–2.4 objects per layout.

## Refits (information only, not used as baselines)
- **Same best subsets, masses refit (non-negative least squares, positions fixed):** k=1 213.0×, k=2 177.1×, k=3 157.1×, k=4 83.4×.
- **Free refit of k clumps** (angles, masses and a shared distance, in the ring plane):
  - k=1: 211.8× (2.7e51 kg at 21.8 Gpc).
  - k=2: 171.4× (2 × 8.6e51 kg at 26.3 Gpc).
  - k=3: 146.0×, but with all negative masses (unphysical).
