# termslice: term-by-term solve of Chris's off-centre CoE geometry per redshift slice

## Method
- Geometry: Chris's compact core Mc plus a uniform shell (Omega_s 0.049, at least 14.3 Gpc on every side of Earth), with Sedov-Taylor flow V = kr (tcore.py wraps ptshell and prof_core.run).
- A = |dV|/separation at each slice's own observation distance Rs (orange 0.21, red 0.47, all 0.31 Gpc), compared linearly with Migkas A/A_all (never H squared).
- Model: A = af(theta; Rs) x [1 - Amp sin(psi + phi)], where theta is measured from the CoE longitude l0.
- Least-squares match to the provided dip+quad terms (C1, S1, C2, S2) for each slice (migkas_zslices_dipquad.json). z<0.03 is ignored as requested.
- D grid: 20, 30, 49.17, 94, 200 Gpc. Both l0 branches were tried, plus a per-cluster distance-distribution variant and a variant solved for the 90-degree-cone-binned profile.
- Mass: Mc = mcc x (4 pi/3) rho_c(H_E) D^3, with H_E scaled so that the curve's mean is <H> = 71.6. H0 = 70 is never an input.
- R2 is scored against the rebuilt Migkas 2020 slice profiles (cycle3/slices2.json).

## Per-slice results (D = 49.17 Gpc)
- Red (z>0.07): Mc/Mcrit 0.00153; l0 15.5 deg; sine Amp -0.031, phase 4.75 rad (dipole peak at l 193 deg); R2 0.74 (the data's own dip+quad fit gives 0.73); Mc 2.1e53 kg.
- All (black): Mc/Mcrit 0.00195; l0 8.3 deg; Amp -0.052, phase 0.29 rad (dipole peak at l 82 deg); R2 0.87; Mc 2.65e53 kg.
- Orange (0.03-0.07): Mc/Mcrit 0.00332; l0 154 deg; Amp -0.140, phase 3.96 rad (dipole peak at l 17 deg); R2 0.93; Mc 4.3e53 kg.
- Solving for the binned profile instead: mcc is about 1.5x larger (red 0.0023, all 0.0029, orange 0.0050), and binned R2 is 0.74 / 0.87 / 0.94.
- The headline Table 1 value (omb, Mc 2.31e53 kg, mcc 0.00171) lies between the red and all solutions.

## Degeneracy
- Every D from 20 to 200 Gpc solves the four terms exactly, so D is not constrained.
- Mc/Mcrit changes by 5% or less across that range, but Mc in kg scales as D^3: from about 1.4e52 kg (D 20) to about 1.3e55 kg (D 200) for red, and 1.7e52 to 1.7e55 kg for all.
- The geometric dipole and octupole are 0.002 or less, so the geometry supplies only the quadrupole and the sine supplies essentially the whole dipole.

## Red vs all consistency
- l0 agrees: 15.5 vs 8.3 deg, a 7-deg difference against about 27-29 deg rms 2l-phase noise (cycle 3c).
- mcc agrees within 1 sigma: 0.00153 vs 0.00195, following quad amplitudes 0.064 +/- 0.031 and 0.084 +/- 0.029.
- The sine dipoles point differently (193 vs 82 deg), but red's dipole is not significant (0.032 +/- 0.031, p_iso 0.83).

## Is the sine plausible as rotation or tilt?
- A solid-body rotation adds only an isotropic term to |dV|/s. It produces no dipole and cannot turn the axis.
- A tilt (peculiar velocity v) gives a dipole of v/(<H> s). That requires about 1030 km/s (red), 1160 km/s (all) and 2060 km/s (orange). These are above the CF4 bulk flow (about 400 km/s) and the CMB dipole (370 km/s), and comparable to Migkas 2021's roughly 900 km/s only for red and all.
- A single tilt must have a fixed direction with dipole falling as 1/s. The slice directions (13, 193, 81 deg) do not fit that.
- A dipole factor multiplying the quadrupole makes 1l and 3l terms. It does not rotate the 2l phase.

## Orange as the black/red model rotated by a free constant dl (one extra parameter)
- Black model: dl = -38 deg (bootstrap median -39 deg; 68% interval +/- 10 deg; distribution heavy-tailed, std 39 deg; 1% of 500 bootstraps give dl >= 0). R2 0.61 rotated vs 0.11 unrotated.
- Black model, binned profile: dl = -40 +/- 10 deg, R2 0.49.
- Red model: best dl = +142.5 deg, which is -37.5 deg for the quadrupole with the dipole flipped (red's dipole peaks at 193 deg, orange's at 13 deg); +/- 11 deg; R2 0.47 (unrotated -0.26).
- This confirms Chris's reading: the orange pattern is the red/black pattern shifted by about -38 deg.
- The remaining misfit after rotation is amplitude: orange's quadrupole is 1.8x black's and its dipole 2.8x.
- The bootstrap +/- 10 deg keeps the model shape fixed. Model-based mocks (cycle 3c) give about 36 deg rms 2l-phase noise in this slice, and 24% of them show a miss this large, so the offset is about 1 sigma against noise.

## What could produce the offset at about 0.2 Gpc
- Axis rotation or tilt (the sine term): cannot produce it. Rotation adds no angular structure, a tilt adds 1l and 3l, and the geometric axis stays fixed to within about 2 deg from 0.02 to 2 Gpc (cycle 3b).
- Local flow: turning the black quadrupole into orange's requires an added local quadrupole with amplitude 0.14 to 0.16 in A and axis at l about 137 to 142 deg.
- That direction is close to the Shapley / Great Attractor / Perseus-Pisces axis (l 307-312 deg, i.e. 127-132 deg mod 180; Perseus-Pisces about 140-160 deg), so a tidal-shear origin is plausible in direction.
- The amplitude is not plausible: it is as large as orange's whole quadrupole, about 5x the dipole fraction a 400 km/s CF4 flow gives at 0.2 Gpc (about 0.03), and it would need coherent density contrast of about 0.3-0.6 on 200 Mpc scales, far above LCDM rms (about 0.05-0.07).
- Most likely: noise plus possibly a modest local shear, rather than a physical rotation of the CoE axis.

## Caveats
- The model is 2D, assigned by longitude only, and uses point directions, while the data are 90-deg cone averages. The binned variant raises mcc by about 1.5x.
- Only the Migkas 2020 L_X-T sample is used. Slices are noisy (see cycle 3c mock p-values).
- The dip+quad target terms come from the provided JSON.
- A is compared linearly. If A scaled as H squared, the implied H anisotropies would halve.
- The mass depends on the <H> = 71.6 convention through H_E. D is unconstrained, so any kg value is conditional on the chosen D.

## Files
- tcore.py, solve.py -> solve.json, solve.log; rotate.py -> rotate.json; make_fig.py -> termslice_summary.png
