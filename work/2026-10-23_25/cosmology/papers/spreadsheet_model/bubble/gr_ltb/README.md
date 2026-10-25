# GR (LTB, dust, Lambda=0) version of the thin shell + core best fit (universe_thin.json "free")
Code reuses papers/cmb_redshift/ltb.py, rays.py, analyze.py (same conventions, G=c=1, Gpc); only the mass function is new (ltb_thin.py).
- setup.py: GR-A = today's areal mass distribution identical to the Newtonian one (point-like core Lc=0.005D, empty interior, uniform
  Omega_s=1.012 shell 79.42-108.02 Gpc, vacuum outside; edges smoothed over +-0.5 Gpc), uniform bang t_B=0, E(r) from R(r,t0)=r, t0 from local H=70.
  GR-B = earlier convention (ltb_cmb case B): same profile shape on the Lagrangian label (scaled by r_obs/D), E=r^2/(2t0^2), t_B=0, r_obs,t0 from R=D, H=70.
- run_gr.py -> gr_results.json, A_GR-A.npy, A_GR-B.npy (A propto D_L^-2 at z=0.05,0.1,0.2,0.3, b=0, CoE at l0+180).
- checks.py -> checks.json (shell crossing, null-condition drift); agelimit.py -> agelimit.json (age limits, bound-Earth window).
Headline: Earth's shell is unbound in GR (2E=+598 GR-A, +343 GR-B); H_r < H_perp at Earth (net mass inside D-R is negative:
core 2.5e54 kg vs hole deficit 5.7e56 kg), so the nearby A(l) peaks across the CoE axis: corr with Migkas -0.71/-0.74,
R^2 (z=0.1) GR-A -67 (no sine) / -66 (sine), GR-B -4.90 / -4.66, vs Newtonian 0.909 (sine), -0.73 (same geometry, no sine), 0.636 (no-sine refit).
Uniform-bang age 13.83 Gyr with this mass today needs 1+2E<0 in the outer shell (2E down to -158): not allowed; max uniform age 11.94 Gyr.
