# CMB shell model — assumptions (cmb_shell.py)
A1 Units c=1, Δt=1. Exactly-c geometry at t=3: observer coasts at v=c, D=v(t+t0) with t0=3 → D=6; R=c·t=3; R_rec=3 fixed = D−R.
A2 Kinematics: everything coasts from a common origin at t=−t0 (v=r/(t+t0)); mass does NOT decelerate anything (gravity enters only as a potential redshift). Inconsistent for large M.
A3 Emission event for direction θ (θ=0 looking at CoE): point at distance R from the observer at t_e=t−R/c=0 (the visible circle). Its radius r_e runs 3 (θ=0) to 9 (θ=180°).
A4 Light leaves at 3000 K "at recombination": the emitting matter crossed R_rec at t_rec=R_rec/v−t0 (from 0 down to −2). The light is taken to cool from t_rec to now by stretch.
A5 Light path: straight line at c through static space, from P(t_e) to the observer's position at reception (D=6). The observer's position at t_e (D=3) is not used for the path. This is a simplification; no expansion is applied to the path itself.
A6 Stretch = cube root of shell volume growth from t_rec to now = (t_now+t0)/(t_rec+t0) (homologous coasting, range 2–6), raised to the power f.
A7 Growth fraction f: share of each velocity that is expansion. Motion through space = (1−f)·v, radial for the source, along the axis for the observer.
A8 Doppler: classical, medium = space frame: ν_o/ν_e = (1 − u_o·k̂)/(1 − u_s·k̂). Undefined when the source catches its own light toward us (denominator ≤ 0) or the observer outruns it (numerator ≤ 0). SR is not used because speeds reach 3c. Undefined directions are excluded from the stats and noted.
A9 Gravity: Newtonian Φ(r) of a spherical density law, with total mass M normalised over its support. Weak-field ν_o/ν_e = 1 + (Φ_e − Φ_o)/c². Flagged as failing when |Φ|/c² > 0.1. Φ at t_e uses homologous scaling Φ_te(r) = λΦ_now(λr), λ=2. Physical units: length unit = 1 Gpc (illustrative choice), so GM/(c²ℓ)=1 at M=2.09e22 M☉.
A10 Density laws: from level_sets_3pct (power law, A+B/r, A+B/r³, quadratic γ=0) on 3–9 only (nothing inside 3 or beyond 9). Step law: flat on 3–9 and 5× that inside R_rec (0–3); the factor 5 is an illustrative choice. Masses 0, 1e18, 1e20, 1e21, 1e22 M☉.
A11 Stats: sky mean weighted by sinθ, dipole = a1/mean from the cosθ projection, residual = rms after removing mean+dipole (reported only; plots show full curves). Azimuthal symmetry is assumed.
A12 Nothing is tuned to 2.725 K, Migkas, or flatness.
