# Volumes and counts: visible universe, our shell, typical void / cluster

Constants are the same as local_web/web_core.py: R_vis = 14.3 Gpc, D = 58.2 Gpc (us to CoE), d = 14.3 Gpc, mean density ρ = 2.9e-27 kg/m³,
1 Gpc = 3.0857e25 m. Packing (PACK) is voids 78% vol / 23.4% mass at δ = −0.70; walls 6% / 8%; filaments 15.5% / 40%; clusters 0.5% / 28.6%.
Real-web cell size is s = 0.05 Gpc (one void per cell, mean cell volume = s³ = 1.25e-4 Gpc³).

## 1. Visible universe
V_vis = (4/3)π R_vis³, R_vis = 14.3 Gpc (the "vision sphere" used in every run) = **1.22e4 Gpc³ (12,249 Gpc³) = 1.2e13 Mpc³**.
(The sphere out to the CoE distance, R = 58.2 Gpc, would be 8.26e5 Gpc³. That is not what we see.)

## 2. Shell we're in (around the CoE, 43.9–72.5 Gpc)
V_shell = (4/3)π[(D+d)³ − (D−d)³] = (4/3)π(72.5³ − 43.9³) = **1.24e6 Gpc³ (1,241,866 Gpc³) = 1.2e15 Mpc³**, which is 101× the visible universe.
The tidal region used in the foam runs (14.3–30 Gpc from us, CoE interior excluded) is 8.82e4 Gpc³.
With the exp-lock density (L = 922.65 Gpc) the shell's mass changes by only −0.24% (1.056e56 kg against 1.058e56 kg if uniform).

## 3. Typical void and cluster
- Void: V_void = 0.78 s³ = 9.75e-5 Gpc³ = **9.75e4 Mpc³**, so R_eff = 28.6 Mpc = 20 h⁻¹ Mpc (h = 0.7). The foam's median was 27 Mpc.
  Mass missing per void = 0.546 ρ s³ = 5.8e45 kg = 2.9e15 M☉.
- Cluster (web junction): there are 24π²/35 ≈ 6.77 junctions per cell (Poisson-Delaunay), so
  V_cluster = 0.005 s³ / 6.77 = 9.2e-8 Gpc³ = **92 Mpc³**, so R ≈ 2.8 Mpc. For comparison, a virial R = 2 Mpc sphere is 34 Mpc³.
  Mass per cluster = 0.286 ρ s³ / 6.77 = 4.5e44 kg = 2.3e14 M☉. The foam's median was 2e14 M☉.

## 4. Counts (N_void = f_void V / V_void = V / s³; N_cluster = 6.77 N_void, which is the same as 0.286 ρV / M_cluster)

|region|volume (Gpc³)|total mass (kg)|voids|clusters|mass missing from voids (kg)|mass in clusters (kg)|
|-|-|-|-|-|-|-|
|visible sphere (14.3 Gpc)|1.22e4|1.04e54|**9.8e7 (~100 million)**|**6.6e8 (~0.7 billion)**|5.7e53|3.0e53|
|tidal region (14.3–30 Gpc)|8.82e4|7.5e54|7.1e8|4.8e9|4.1e54|2.2e54|
|CoE shell (43.9–72.5 Gpc)|1.24e6|1.06e56|**9.9e9 (~10 billion)**|**6.7e10 (~67 billion)**|5.8e55|3.0e55|
|sphere R = 58.2 Gpc (for reference)|8.26e5|7.0e55|6.6e9|4.5e10|3.8e55|2.0e55|

Net mass is zero. The void deficit (0.546 M) equals the added mass in clusters (0.286 M), filaments (0.19 M) and walls (0.07 M).

## How this compares with the earlier "~10⁹ voids, several billion clusters"
That estimate was for the **tidal region (14.3–30 Gpc)**, the 7.1e8 in summary_table.md. The full CoE shell is 14× bigger,
so it holds about 1e10 voids and 7e10 clusters.

## Caveat on clusters
The foam puts all 28.6% of the cluster mass into knots averaging 2.3e14 M☉. Counting only observed rich clusters (≳1e14 M☉, roughly 1–2e-5 per Mpc³)
gives about 2–5× fewer: ~1–2e8 in the visible sphere and ~1–2e10 in the shell.
