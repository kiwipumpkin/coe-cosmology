# Phase 5: picking d with the CMB TT chi2 (IN PROGRESS)
Held fixed: <H0> = 68.64 km/s/Mpc and age 13.787 Gyr (5% penalty). Migkas A = af^2/<af^2>.
- fitd5.py d: refits D, mc, l0, eta and k at fixed d with unbounded Nelder-Mead on log parameters. Only the physical guards in k4b.sky apply. Checkpoints: fitd_<d>_ckpt.json. Launched by run_fits.sh (nohup).
- chain5.sh (nohup): waits for all 9 fitd_<d>.json files, then runs cmb_one.sh once for each unique D (rounded to 0.05 Gpc), then runs summary5.py and the plots.
- p5_comp5.py --obs 1: compensation targets the observed Planck low-l sky (l=2..20 at full weight, cos^2 taper to l=40) instead of zero, so C2 is preserved.
- p8_scratch5.py: Dfix used for l=2..128, not 2..40. The l=40/41 truncation step caused the l~30-45 spikes.
KEY CAVEAT: the CMB chain (ltb5, p2_bg5, p5, p8, p6) depends only on D, not on d, eta, k or mc. So TT chi2 "vs d" is really TT chi2 vs the refit D(d). D(d) spans only 27.33-27.44 Gpc.
Interim Migkas results (best start per d, read from the ckpts): chi2 8.516-8.521 for every d, so Migkas is flat in d. log10 k runs away to 1e49-1e244, so k is unconstrained and runs off to infinity.

## Phase 5 (NEW, 02:30 CT Oct 26): CMB built from the Migkas density profile
- bg_mig.py: the background is built from rho(r) = exp(-(r/L)^eta) (Migkas profile; d and k enter through L via the CatWISE lock), with <H0> = 68.64. frac is SOLVED so that H t0 gives age 13.787 Gyr; frac/mc is reported as the cross-code check. Normalisation is documented in the script header.
- Validation (poly mode, SURF=rho), bg_VALpoly.json: age 13.7866 Gyr, R_ls 21.7395 Gpc, matching phase4b (21.74). That covers the background only; the full TT validation (VALpolyA, SURF=a) is still running.
- famP_mig.py: famP with the Migkas profile option. fitdk5.py: Migkas fits at fixed (d, log10 k) = {0.5, 2, 6} x {-1, 0, 1}, still running.
- chain6.sh (nohup, OMP threads = 1, launched 02:25 CT): queue6.py -> cmb_mig.sh, 3 runs in parallel. After all 9 fitdk fits finish it re-queues the J_* tags, then runs summary6.py and the plots (chi2_vs_d6, TT_best6 linear/log, each with its .py). Resume-safe: re-run `nohup bash chain6.sh &`.
- The old D-only rebuilds (chain5) are NOT run.

## Phase 5 chain6 failure diagnosis (NEW, Oct 26 ~13:30 CT)
- chain6 and its p5_comp5 jobs (K1, K3, A_d0.5) stopped. Invalid outputs for A_d0.001, A_d10, K-1 (p5_bg, p8, bg map/trace/json) moved to bad/.
- The bg maps were already finite; the NaN/inf is downstream. Root cause is the background itself, not a numerics bug:
  the Migkas-profile backgrounds give an O(1) CMB anisotropy. Example K-1 (rebuilt): observed 1+z ranges 1489-6403 over the sky,
  map rms 0.46, rms after removing the dipole 0.35, l=2 amplitude 0.13. Phase4b/poly: rms after removing the dipole 8.6e-4 (VALpolyA 2.8e-3).
  Compensation with MMAX=1e22 then leaves a residual of ~0.08-0.17 rms (27-38% of masses at the bound), Dfix ~1e11 uK^2, so p8 goes inf/NaN and the comp lstsq SVD fails.
- Bug fixed on the way (famP_mig.py): te(x) was tabulated only on x in [0.6,1.6] and linearly extrapolated, but rays hit at x~0.2-2.6.
  The table is now widened until it covers every hit (cubic, no extrapolation), and the script aborts if any ray misses. This did NOT remove the O(1) anisotropy.
  For K-1 the hit surface jumps from x~1.85 to 0.67 between psi=110 and 115 deg.
- cmb_mig.sh now fails loudly: every step is checked with check_finite.py and prints ONE_FAIL instead of continuing to ONE_DONE. chain6.sh flags any ONE_FAIL.
  Backups: cmb_mig_pre_fix.sh.bak, chain6_pre_fix.sh.bak, famP_mig_pre_fix.py.bak.
- Chain NOT relaunched: every Migkas-profile bg map (all A_*, K*, Tk*) has the same O(1) amplitude (min ~-0.8, max ~+2), so it would produce no valid TT.

## Phase 5: recombination-surface definitions tested (NEW, Oct 26 PM CT)
- famP_mig.py: added SURF=t (constant proper time since the common bang, t_rec = time our shell had a=1/1100; override with env TREC), and PSI env (custom ray list). Existing: SURF=a (each shell's own a=1/1100), SURF=rho (density = our shell's density at a=1/1100). Backup: famP_mig_pre_surf.py.bak.
- Runs in surf/ (run_surf.sh, metrics.py -> metrics.json, fine.sh fine psi 100-120 scan, rho_diag.py -> rho_diag.txt). New bg_*_St / *_Srho files in this dir.
- Results (dipole-removed rms of the bg map; l=2 rms):
  K-1  SURF=a 1+z 1489-6403, rms 0.46, nodip 0.35, l2 0.13, discontinuous(psi 110-115)
  K-1  SURF=t 1+z 1008-30840, rms 1.18, nodip 0.70, l2 0.33, discontinuous(psi 113-115)
  K-1  SURF=rho FAILS: no shell-wide time where rho equals our shell's value (brentq sign error) -- surface undefined
  A_d2 SURF=a 1+z 1262-21558, nodip 0.57; SURF=t 1003-30650, nodip 0.70; SURF=rho FAILS
  VALpoly(A) SURF=t 1092-1106, nodip 8.7e-4, l2 5.5e-5, continuous (SURF=rho 8.6e-4, SURF=a 2.8e-3). R_ls 21.7395 / 21.7439 Gpc unchanged by SURF.
- Root cause (rho_diag.txt): with common bang, the Migkas shape is the INITIAL mass distribution, so at t_rec density is already O(1e3..1e-95) inhomogeneous:
  K-1 rho(t_rec,x)/rho(t_rec,1) = 1950 at x<=0.5, 3.9e-8 at x=1.3, 1.8e-95 at x=2.5 (poly: within 0.1% for x<=2). Local log-slope at our shell 43 (K-1), 1.2e4 (A_d2).
  Rays leaving outward cross near-vacuum (no plasma: tau never reaches 1 there), rays inward enter the dense core; the 'jump' is the steep but continuous sweep of rays across the core edge near psi~113-115.
  => No recombination-surface definition fixes this; the background itself violates the LTB-literature assumption of a near-homogeneous early universe (no decaying mode). Full comp/TT run (step 4) NOT done.

## Phase 5 routes A/B (NEW, Oct 26 ~15:00 CT) -- routes/
- Route A (fitroute.py A EPS): common bang, growing mode, Migkas profile = PRESENT-DAY density (homog.py, dist2 physics); hinge penalty forces max|rho(t_rec)/rho(t_rec,D)-1| <= EPS over x reached by rays. Route B: local log-slope at x=1 locked to CatWISE 0.04104 (pins L/D given eta).
- Results (metrics: routes/metrics_routes.json; all age 13.787, R_ls 21.744-21.745 Gpc, no shell crossing):
  A eps=0.1 : chi2 8.706 r2 0.752 D 27.88 mc 0.137 l0 19.06 eta 0.517 L 1.67; implied local k 54. St: 1+z 1051-9863 nodip 0.167 l2 0.080 discont.; Sa: nodip 0.258
  A eps=0.01: chi2 10.14 r2 0.730 D 30.42 mc 0.841 l0 20.19 eta 0.0946 L 3.0e-5 (LD at 1e-6 = FLAGGED, numerical floor); fit age 12.43 (penalty). St: 1+z 1091-1292 nodip 0.0215 l2 0.019 cont.; Sa: nodip 0.0374
  B         : chi2 18.10 r2 0.471 D 25.54 mc 0.421 l0 28.00 eta 1000 (FLAGGED at k4b guard 1e3 -> step profile). St: 1+z 1100-29869 nodip 0.634; Sa: nodip 0.286
- None reach nodip ~1e-3 => full comp+TT+p8 NOT run. A eps=1e-3 fit running (nohup, PID 1795734); routes/chainA3.sh (nohup) then builds its bg + metrics.
- density_routes.png (+ density_routes_plot.py). v1_linM/: superseded first attempt (linear M interp, ~3e-3 spurious delta).

## Phase 5 web test (NEW, Oct 26 ~15:30 CT) -- web/
- web/run_web.sh: p5_comp5.py (cmb_mig.sh settings: NB 3000, LC 40, LP 60, --obs 1 Planck low-l target) on bg_RA_{0.01,0.1}_St, MMAX = 1e22 (OUR bound) / 1e24 / 1e27. Outputs p5_web_RA<E>_M<M>.{json,npz} in phase5/.
  (First launch ran 6 in parallel; the eps=0.01 jobs were OOM-killed and rerun sequentially.)
- web/dTT_rms_vs_angle_web.py -> web/dTT_rms_vs_angle_web.png + web/web_metrics.json. delta = M/(rho_crit(68.64) * shell volume/NB).
- eps=0.01: bg rms 0.0509 / dipole removed 0.0215 / l2 0.0194 -> cap 1e22: 3.02e-5 / 3.02e-5 / 3.5e-6; comp max 3.7e-4 rms 2.6e-5; 0% at cap; |delta| max 0.89 (no void below -1). PHYSICAL.
  cap 1e24: rms 1.64e-5, comp max 1.0e-5 rms 8.7e-7, |delta| max 1.22 (5 voids < -1). cap 1e27: rms 1.64e-5, |delta| up to 54 (673 voids < -1) = unphysical.
- eps=0.1: bg nodip 0.167 -> cap 1e22 2.97e-3 (|delta| 3.2, 92 voids < -1); 1e24 7.4e-4 (|delta| 40); 1e27 2.3e-4 (|delta| 1e4). Not physical below ~3e-3.
- Caveat: bodies fit l<=40 (penalty to 60) on Nside 64 maps; rms is at Nside 64 (l<=191).
- web/tt_web.sh (nohup, resume-safe, fail-loud): p8 -> p5_fld -> p6 -> tt_web_summary.py (web/tt_web.json) -> TT_web_plot.py / TT_web_plot_log.py for web_RA0.01_M1e22. Launched ~15:35 CT.

## Phase 5 web seed robustness (web/run_seeds.sh, web/seed_rms.py)
Route A eps=0.01, cap 1e22, settings identical to run_web.sh, --seed varied (11,23,101,777,1234,3141,5000,9001 + original 2000; n=9).
Outputs p5_web_RA0.01_M1e22_s<seed>.{json,npz}, web/seed_metrics.json, web/seed_rms.png.
dipole-removed rms: mean 2.87e-5, std 2.1e-6, min 2.59e-5, max 3.12e-5. l=2 amp 3.45e-6 (std 1e-8). comp max 1.6e-4..4.3e-4, comp rms 2.06e-5..2.67e-5.
frac at cap 0 for all seeds; max|delta| 0.55..0.89; delta<-1 count 0 for all seeds. Cancellation is robust to body placement.

## Web TT bump ablation (web/ablation/, Oct 26 2026)
Baseline web_RA0.01_M1e22. Model TT = Df (p8 stochastic field) + Dfix (l=2..128 = Dl of bg+bodies residual map). `ablation.py` is a frozen-field
ablation: it uses the baseline p8 best-fit params and CORR without refitting, recomputes Dfix per variant (it reproduces the baseline χ² 27001.57 exactly), and fails loudly on NaN.
ISW off = acoustic_model.EX.KI zeroed in memory (no files edited). "Near" = the nearest 10% of bodies from the observer (300 bodies, r ≤ 9.10 Gpc;
all have a softening core >5°). Outputs: ablation_metrics.json, ablation_spectra.npz, TT_ablation_log.png/.py, diag_bins.py/.json.
Result: no ablation removes the bump. No web removes it at l≥15 but leaves the raw background (C2 6.7e9). ISW off and nearest-10% removal leave it unchanged.
diag_bins: the l=15-35 power comes from the outer bodies (0.75-0.875 R, next to the CMB shell). Removing them takes D_24 from 9710 to 497,
while bodies within 0.5 R contribute ~0. The bump is the bodies' own small-scale imprint, i.e. a compensation error (Resid−Tobs: D_24 = 8641).

## Inverse problem: uniform random web + minimal displacements (web/inverse/, seed 1, Oct 26 2026)
Setup (inverse2.py): 3000 bodies uniform in volume in [0.01, 0.975] R (R = 21.7439 Gpc), lognormal delta (sig 0.078, delta in [-0.220, 0.318]),
1406 clusters / 1594 voids, masses fixed, on Route A eps=0.01 St background. Greedy moves (25/batch) to random candidate positions + LM refinement
of moved positions; objective = weighted alm residual vs Planck SMICA low-l target (l=1..40 fit, 41..60 penalised). MAXMV=1000 is a compute budget.
Run: `python3 inverse2.py --seed 1 --MAXMV 1000` (via run_inv.sh), stopped on STALL at 938 moves (869 unique bodies = 29.0%; 26.6% of |M|).
Baseline: the pure random web does NOT cancel anything (as expected): dipole-removed sky rms 2.36e-2 vs background alone 2.15e-2.
Result: dipole-removed residual rms (sky − target) 2.36e-2 → 9.8e-3 (100 moves) → 8.8e-4 (400) → 4.70e-4 (final, plateau from ~800 moves).
Earlier shell-concentrated free-mass web reached 2.55e-5, i.e. this inverse solution is ~18x worse and ~29x the Planck target rms (1.64e-5).
cancel max |resid| 7.7e-3. Displacements: median 8.2 Gpc, mean 12.0, p10 1.6, p90 28.2, max 42.1 Gpc.
Radial: moved bodies end at mean r/R 0.904 (originated at 0.807); 93.0% end beyond 0.75 R vs 54.5% for uniform -> they DO pile up near the shell.
132 moved bodies sit at the r = 0.975 R geometric limit (RMAXF, one softening length inside the shell) — FLAG: pinned at a bound; 27.5% within 1% of it.
TT (tt_inverse.py, FROZEN field = baseline p8 params, not refitted; l=2..128 replaced by D_l of the bg+web sky): chi2 total 1.734e9, l2-30 1.161e7,
l30-2500 1.725e9; r2 total 0.0034. Residual sky power is huge at l~20-130 (D_20 7.9e4, D_50 1.45e6, D_100 3.6e5 uK^2 vs Planck <~6e3): the l~15-35 bump
is present and much worse, extending to l~130 (the l>60 part is outside the fitted range and uncontrolled). C2 = 167 uK^2 (Planck 226).
Peaks: find_peaks returns 103-115 (noise in the sky term), not acoustic peaks — the l<=128 band is dominated by the web residual; above l=128 the frozen field is unchanged.
Files: inv_s1.json/.npz, tt_inverse.json/.npz/.log, TT_inverse_log.png/.py, radial_moved.png/.py, moved_map.png/.py (clusters vs voids, from->to, mass),
rms_vs_moves.png/.py/.json (dipole-removed residual from replay snapshots + log cancel_rms), rms_vs_angle_stages.png/.py/.json (0/100/400/final).
inverse2_replay.py = identical deterministic rerun saving per-cycle snapshots to replay/ (reproduces inv_s1.log cycles to ~1e-13); rerun rms_vs_moves.py when it finishes.
Seeds 2,3 continue sequentially under run_inv.sh (nohup, per-seed checkpoints).

## Inverse problem v3: mass + position free, fit to l<=130 (web/inverse/v3/, Oct 26 2026) -- QUEUED, NO RESULTS YET
- `v3/make_target.py` -> `target_Tobs_l130.npy`: SMICA alm (native Nside 2048) l=2..130, synthesised at Nside-64 pixel centres
  (no ud_grade pixel-window loss; check: Dl(target) = SMICA Dl to <0.5% at l<=128, e.g. l=20: 720 uK^2, l=100: 2746). rms 2.56e-5
  (the v2 l<=40 target had rms 1.64e-5, so "Planck level" differs between v2 and v3).
- `v3/inverse3.py` (new file; ../inverse2.py untouched): same uniform random 3000-body start as v2 for a given seed (identical rng order).
  Objective: alm residual l=1..130 vs target. Each body's delta free with position; only physical bound delta>=-1; NO upper cap.
  Geometric limits r in [0.01,0.975]R (0.975R = one softening length inside the shell). Greedy L0 on bodies changed (move+optimal
  mass, or mass-only), LM on (position, delta) of the changed set every 25 changes. Pinned counts (0.975R, RMIN, delta=-1) in JSON.
  Per-cycle checkpoint (small; matrices rebuilt on resume) + snaps/. Fails on NaN. Smoke test (NB=200, NC=500; not a result) in v3/test/:
  resume and exact incremental-map drift ~1e-16 verified.
- `v3/tt_inverse3.py` (frozen field, as ../tt_inverse.py, adds D_l at 20/50/100), plots `rms_vs_moves.py`, `radial.py`, `moved_map.py`,
  `TT_log.py` (each PNG has its .py beside it).
- `v3/run_v3.sh` (nohup, log run_v3.log): waits for the v2 chain (run_inv.sh) to exit, then seeds 1,2,3 one at a time, TT + plots after each.
