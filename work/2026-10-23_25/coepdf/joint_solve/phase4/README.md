# phase4 (NEW, 2026-10-25 CDT): <H0> locked, age anchor, gradient free
Dynamics code: prof_core ('eta' bubble) for the Migkas sky AND the age (age = Ht*977.8/68.64 Gyr).
famP/dist2 LTB (CMB background) is matched to it by equal H*t0: frac=0.0185 gives Ht0=0.9679 (13.79 Gyr) -> R_ls = 0.6105*c-units = 21.77 Gpc (ltb_scan.py, ltb_f*.log).
Files:
- pc_probe.py/json: prof_core Ht vs mc (age 13.8 Gyr needs mc of about 0.06-0.2).
- k_scan.py: model functions (its grid run was aborted, k_scan_aborted.log; superseded by the joint fit).
- joint_fit.py -> joint_T0.05.json: differential evolution over logD, logd, logmc, l0, logeta, logk (k = delta/0.031), age penalty ((age/13.787-1)/0.05)^2.
- profile_k.py -> profile_k.json; profile_plot.py/.png/.json: profile chi2 vs k, with (A) d free and (B) d = 14.30 Gpc.
- p2_bg_build.py (D H0 FRAC) -> p2_sm_D27.92_H68.64_f0.0185_ns64.npy
- p5_compensate.py -> p5_a_age*.json (MMAX scan 1e14/1e18/1e22, LC 24)
- p8_scratch.py: stochastic field from scratch (random draws in the bounds, no earlier parameters) -> p8_age.*
- p6_tt_compensated.py -> p6_age.*

## Results (finished 2026-10-25 evening CDT)
- joint fit: D=27.92 Gpc, d=2.00 Gpc (AT LOWER BOUND), mc=0.0651, l0=18.81 deg, eta=2.20, L=10.8, k=99.9 (delta=3.10), age 13.787 Gyr, R_ls=21.77 Gpc. Migkas (A=(H/<H>)^2, n=71): chi2 8.57, r2 0.757; flat 45.14.
- profile k (A, d free): chi2 is monotone decreasing in k, min 8.546 at grid edge k=1000; 1sigma lower bound k>2.84 (2sigma k>1.03); k=1 dchi2=+4.03. No upper bound. (B, d=14.30): k>16.8 (1sig), k>6.8 (2sig); k=1 dchi2=+35.
- p8 (16 starts, converged/finished): TT chi2 l=2-2500 ensemble 12470.9 (l2-30 64.5, l30-2500 12407.4); seed42 13037.4. LCDM 2558. Peaks l=246,534,828,1124,1422.
- compensation with field (p5_age_fld, NB3000, MMAX1e22): max|dT/T|=2.99e-5, rms 1.03e-5 -> FAILS 1e-5. NB8000: max 2.58e-5, rms 9.4e-6 (residual dominated by l=2 ~9.1e-6) -> FAILS. RMAXF0.95/NB8000 run OOM-killed. Background-only best: 2.4e-5 (M1e22).
- p6 compensated seed42: l2-2500 13028.6.
