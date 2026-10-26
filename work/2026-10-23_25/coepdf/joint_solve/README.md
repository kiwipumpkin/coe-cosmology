Joint solve (2026-10-25). Run order: s1_migkas_baseline.py; s2_shear_compensate.py <seed> <N> <Lc> <Mmax> <Rmax> <lam>;
s3_tt_refit.py <s2 npz> <maxfev>; s4_lowl_fraction.py <s3 npz>; s5_migkas_check_and_figs.py <s2 npz> <s3 npz> <s4 json>.
Best set: s2_s2000_N3000_L12_M3e+16_r3.5_lam1e-05 (+ matching s3_/s4_ files, s5_joint_summary.png, s5_migkas_check.json).
s3_..._M1e+17_r3.0_lam1e-03.* : TT refit of an earlier lam=1e-3 compensation whose s2 npz was overwritten (tag lacked lam); kept for record only.

## Phase 2 (phase2/) — consistent background, TT extended search, actual-body l=2 cancellation
- p1_migkas_Dprofile.py: Migkas refit vs D (30/34.96/42/50/58.2): chi2 9.528/9.250/9.594/10.277/10.886 -> D=34.96 chosen (eta at 0.03 edge; NM loose).
- p2_bg_build.py: background map from famP with D as parameter (D=42 reproduces hl_R4 sm, corr 0.99999999); map ~independent of D.
- p3_background_relations.py: D=34.96 -> LTB H0=34.31 km/s/Mpc, age 0.314 Gyr, R_ls=0.18957 Gpc.
- p4_tt_extended.py (15-param P(k), 15 starts x 3000 fev; still running at writing) ; p4b_eval_ckpt.py: interim best of finished starts -> p4b_interim.json/npz.
- p5_compensate.py: bodies at R=0.18957 (a), R=4 legacy (b), seed-42 field l=2 cancel f=0.26 with s3 field (c) and p4b field (c_p4b). No body at bound.
- p6_tt_compensated.py: TT chi2 with compensated seed-42 low-l sky (l<=40 from map). p6_s3.json, p6_p4b.json.
- p7_migkas_bodies.py: Migkas check with bodies -> p7_migkas_bodies.json. p7_summary_fig.py -> p7_summary.png.

## rev7 additions
- phase2/p8: final field fit; p5_c_p8: compensating bodies for p8 (final); p6_p8: TT/low-ell; p7: Migkas re-check with bodies (key p5_c_p8.npz); p5_err: true-error check.
- fig7/N4_tt_stats.py -> N4_tt_stats.json (LCDM comparison, peaks); fig7/N5_seed_distribution.py -> N5 json/png (1000 realisations, rng 2026).
- Fig04_amp_full.py now reads p7 key p5_c_p8.npz (old key p5 run a no longer in p7 json).
