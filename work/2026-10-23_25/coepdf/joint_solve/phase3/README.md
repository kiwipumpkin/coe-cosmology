# phase3 — H0 fixed to observation, sqrt(A) relation, full-model refit  (NEW, 2026-10-25, CDT)

## H0 (h0_sources.py -> h0_sources.json)
Planck18 67.4±0.5; DESI DR2 68.51±0.58; CCHP TRGB 70.39±1.94; SH0ES 73.04±1.04; TDCOSMO 71.6+3.9/-3.3;
masers 73.9±3.0; GW170817 70+12/-8 (ACT DR6 68.22±0.36 listed, excluded: shares Planck data).
Inverse-variance mean 68.64±0.35 (scatter chi2 28.7/7). Band = span of ±1σ for σ≤2 sources: 66.90–74.08.

## H0 in the model
dist2 LTB is self-similar in u at fixed frac: H0 = 1.00016·u·(c/D), c/D=8.5753 at D=34.96 (ltb_scan.py).
H0 is imposed exactly via env U; Migkas sky shape is u-independent, so chi2/r2 identical at 66.90/68.64/74.08.

## sqrt(A) (audit E13) — the CSV issue
catwise-eta-migkas-fit-d14.3.csv column A_over_Aall is the H-like a_f=|v_rel|/r (prof_core l.75) fit as A∝H.
migkas_refit_sq.py refits {D,mc,l0,eta} with A=a_f^2/<a_f^2>: D=34.748, mc=3.876, l0=20.14, eta=0.0300,
chi2_A=9.478 (flat 45.14), r2=0.741, n=71. Curve: migkas_refit_sq_curve.csv.
eta_probe.py: eta<0.03 has no L root under the delta=0.031 lock -> bound is physical.

## N6 (NEW, layout unchanged): N6_H0_direction.py/.png/.json, N6_H0_r2.py/.json — curve = refit H/<H>, no extra sqrt.

## CMB chain at each H0: p2_bg_build.py (bg map), p5_compensate.py (bodies), p8_field_H.py (field), p6_tt_compensated.py (TT chi2 vs LCDM 2558).

## Age (three numbers, three different things)
famP/dist2 LTB: t0 = Ht0/H0, Ht0=0.0110 at frac=20000 -> 0.157 Gyr at 68.64 (set by frac only).
prof_core: age = Ht × (1/70 km/s/Mpc Hubble time, 13.97 Gyr; hard-coded H0=70 in bubble_core.TH0);
  Ht=0.199 (old locked params) -> 2.78 Gyr; Ht=0.349 (squared refit, mc=3.876) -> 4.88 Gyr (4.98 if rescaled to 68.64).
prof_core and famP are separate dynamical models (2-D bubble with mass parameter mc vs LTB with frac); they are not tied together.

## Results at H0=68.64 (R_ls=0.09475 Gpc)
p8 field (12 starts + 3 polishes): TT chi2 l2-2500 = 14521 (mean), 15166 (seed42), 15169 (seed42 compensated) vs LCDM 2558;
phase2 at R=0.18957 gave 2616. The field fit does NOT recover at the halved R_ls -> main thing that breaks.
p5 with field (f=0.26): max|dT/T| resid 9.95e-6 (<1e-5, marginal), rms 4.64e-6, D2 sky 45.1 = target.
