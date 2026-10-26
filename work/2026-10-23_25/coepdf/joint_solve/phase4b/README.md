# phase4b (NEW, 2026-10-25 CDT): no arbitrary d floor. <H0>=68.64, age 13.787 Gyr (5% penalty), k free.
Bounds removed/widened. Old (phase4): joint_fit d>=2 Gpc (arbitrary), logk<=2 (k=99.9 sat at it), D<=80; k_scan.sky guards 1<d<60, eta<8, mc<80, L/D<3.
New (joint_fit4b.py, k4b.py): D 2-316, d 1e-3-60, mc 1e-5-1e3, l0 0-360, eta 1e-4-1e3, k 1e-3-1e5; physical guards only: d<D/1.05, L root exists, vrmin>0.
(At d<~0.1 Gpc with eta~2, k~100 the shell collapses (vrmin<0) -> rejected physically; other eta allow small d.)
## Results
- joint4b.json (DE 120 gen + NM): D=27.34, d=0.001 (AT NEW FLOOR; chi2 flat in d: 8.565 @1.33, 8.54 @0.43, 8.52 @<0.03), mc=0.0560, l0=18.64, eta=11.5, L=12.42, k=1260 (delta=39), age 13.79 Gyr. chi2 8.517, r2 0.759 (n=71, flat 45.14).
- profile_k4b (A, d free): chi2 8.549 (k=0.1) -> 8.517 (k>=100): total range 0.032, k UNCONSTRAINED (dchi2<<1 over k=0.1..1e5). eta absorbs k (3.1->12.4), D 27.74->27.33.
  (B, d=14.30): incomplete at report time (k=0.18: 41.1, k=0.32: 31.4); phase4 B gave k>16.8 (1sig).
- ltb_scan4b (DSHAPE=27.34, frac 0.0185): Ht0=0.9678, eta_rec_D=0.7953 -> R_ls=21.74 Gpc.
- p2 background: p2_sm_D27.3385_H68.64_f0.0185_ns64.npy (rms 2.42e-3).
- chain4b.sh (running at report): p5 bg-only (LC24, WH1) -> p8 from scratch -> p5 with field, WQ=1,5 (quadrupole-weighted; memory-light) -> p6 TT.
