# Feedback on the four PDFs (25 Oct 2026)

## 1. Rev6 workflow paper (52 pp)
- Reads as an honest methods notebook. It already adopts the audit's E13 fix (§4.1, p.12: H/H_ref = √(A/A_all)) and reports the fitted Fig 4 dipole as two parameters, with R² before and after (A 0.446→0.712, B 0.483→0.785).
- Still inconsistent: D = 58.2 Gpc in the abstract, §2.2 and §5.1; the free η-edge best fit at D = 34.96; and D = 42/39.1 in the §5.6–5.9 CMB runs (Table 7, p.34). Rev7 cleans up most of this, but main_rev7.tex still has 6 mentions of "58.2" and 7 of "42 Gpc", so each needs labelling as historical or removing.
- §2.2 admits the particle horizon (≈12.7 Gpc) is smaller than d = 14.3 Gpc. This is still unresolved.

## 2. Referee report on rev5
- Fair and well calibrated. The key point is §2a: the A-curve shape, the rings and the gradient sign were predictions made before the data. That is the strongest card, but it needs **dated evidence** (preprint, notes, commits).
- Right: the χ² over correlated bins (M1); a parameter count against dipole+quadrupole (M2); the CatWISE lock ignoring the kinematic dipole (≈0.007), bias and magnification (M3); the moving-lens frame (M8).
- Points where it now overreaches or is out of date: M6 (A_F Fourier fit) was withdrawn by the reviewer. The B "near-null" claim was corrected.
- On the CMB it treats the work as an existence demonstration, which matches your framing.

## 3. Math/code audit
- Part A checks most of the algebra as correct. The real errors are small: E1 (2.06e-9, not 2.2e-9) and E4 (L = 920.3, not 922.6).
- Structural: E2 (one background; flat Ω_m=1 with 13.8 Gyr gives H₀ = 47.2, which clashes with the H₀ = 70/Ω_m 0.3/71.6 values used elsewhere) and E6 (42 Gpc hard-coded in the CMB code against 58.2 in the cluster fits).
- E8: the l₀ values are 15.2/15.6/19.9/100/0°, and there is a CoE/anti-CoE sign flip between redweb2 and physics.py.
- E13 (A ∝ H², so H ∝ √A) is physically right: Migkas constrain A·H₀².
- **Spot check on our N6:** N6_H0_r2.py does take √ of the curve in catwise-eta-migkas-fit-d14.3.csv. The H-map report traces those law curves to a_f = |v_rel|/d, which is an H-like quantity that the old code compared directly to A. If that CSV column is the H-like a_f and not a model A, then our √ halves the amplitude a second time. In that case the 9.8% peak-to-peak and r² 0.744 would be wrong. This **must be checked before rev7 quotes N6 numbers**.
- B1/B2 (shared bulk velocity should cancel; a hidden 100% P₂ power modulation) are the most serious findings. Rev7 §open already lists them as unresolved, priority 1.
- B5: the hand-inserted ℓ=220 hump in the old Fig 20 must be captioned. B9: the seed-42 single-draw statistics.

## 4. 2-D H sky map vs Migkas
**What it finds:** the full-sky axis fit at fixed amplitude (s=1) puts the anti-CoE at (11°,+7°). That is 8° from our 1-D axis (15.2°,0), with Δχ²=0.6 to move it, and the Δχ²≤2.3 radius is ≈14°. The sense matches physics.py; the redweb2 sense is strongly disfavoured (Δχ² ≈ +14–18).
- The CoE law at fixed amplitude has the best ΔAIC (−2.5) and ΔBIC (≈−14 against dipole+quadrupole).
- But the shuffle null with a full-sky axis search gives p≈0.4. Even the anisotropy itself is typical (p=0.61 at 41 modes; p=0.036 only in the leading 10 modes).
- The free amplitude s=0.92±0.38 cannot tell s=1 from s=½. The best free fit is a sign-inverted "cold pole" at l≈280°, only Δχ²≈2 better.

**Methodology:**
- Sound points:
  - The shuffle keeps positions, T and errors, and includes look-elsewhere.
  - The forward model is checked by injection.
  - The null is calibrated (median χ²_iso ≈ k).
- Weak points:
  - The 75° cones leave only ≈4–6 effective modes (N_eff), so the test has little power by construction. A null result is expected even if the effect is real.
  - The 99%-variance truncation (k=41) dilutes signal. The result moves with k, and the OLS axis wanders at low k.
  - Its estimator does not reproduce the Migkas maximum (1.056 at (17,−9) against 1.074; its global maximum is at the l≈180° limb).
  - The law curves are single-depth (d=14.3 Gpc), not averaged over cluster depth (z<0.3).
  - No catalogue-rotation null was run.
  - The author admits the effective parameter count is >2 because the laws were tuned on the same clusters, so the AIC edge is optimistic.
- It also notes that Migkas's own 3.6σ has no look-elsewhere correction. That point is fair and cuts against the data's anisotropy, not only against CoE.

**For your framing:** this fits a "toy universe, right direction" paper well. The axis, the sense and the economy all point the right way, and nothing is significant. Say exactly that. Don't claim detection. The (11°,+7°) result is independent support for the axis, and for the physics.py sign over redweb2.

## Conflicts between documents
- D: 58.2 (rev6, referee, audit), 34.96 (rev6 best/rev7), 42/39.1 (CMB code). The H₀-anchored refit gives D ≈ 17.1, which breaks the CatWISE-locked Migkas fit (χ² 17.67).
- The CoE/anti-CoE sign: the H-map work settles it in favour of physics.py.
- Amplitude: the audit says the model A anisotropy is "half what it should be", but the H-map's free s does not need s=1.

## Priority actions for rev7
1. Verify whether the law CSVs are H-like or A-like, and fix the √ in N6 and the H₀ refit if needed.
2. Use one D, one background and one axis convention (physics.py sense) throughout. Label or remove the residual 58.2 and 42 values.
3. Add the 2-D map result as a short section: the axis (11°,+7°), sense, ΔAIC/BIC, and p≈0.4 stated plainly as qualitative.
4. Keep B1/B2 as the headline open CMB problem. Run the relative-velocity kernel and compute the P₂ modulation.
5. Add a dated predicted-vs-fitted table (referee §2a).
6. Small fixes: E1, E4, E9 (call the contours slices), E10 (one peak finder), and the B5 caption.
7. Optional: depth-averaged laws and a rotation null in the 2-D test.
