# Heavy-core / R = Migkas-distance refit in exact LTB (ballpark)
Model A = Chris's definition: directional expansion rate (|dV|/separation analogue = H(psi) at fixed z, i.e. D_L^-1, LINEAR, not H^2),
divided by A_all = mean over the 360 sampled 1-deg directions; compared with Migkas A/A_all. Headline R^2 = with the sine (Amp, phi).
- rm_core.py: evaluator (layout: core Lc=0.005D + empty interior + uniform shell D-R..D+R + vacuum; GR-A, t_B=0, H_loc=70), l0/Amp/phi fit.
- grid_mc.json / flipsearch.json: core-mass scan at R=14.3 Gpc, Omega_s=1.012 (NOTE: run with A propto H^2 before the correction).
- rgrid.py -> rgrid.json: R = 0.316 Gpc (median z=0.075) and 1.19 Gpc (z=0.3), grid over D/R, Omega_s, Mc/Mcrit (A linear).
- plot_best.py -> gr_fig2_R_migkas.png / _web.jpg, best_rows.json.  iso_newton.py -> iso_newton.json (isolation tests).
