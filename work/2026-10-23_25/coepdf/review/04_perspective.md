# Peer Review Report — Reviewer 3 (Perspective)

## Reviewer Information
- **Role:** Perspective (card #4). **Identity:** inverse-problems / applied mathematician. **Focus:** identifiability, model flexibility, what the fits can and cannot show.

## Overall Assessment
- **Recommendation:** Major Revision. **Confidence:** 3/5. **Calibration:** NOT_CALIBRATED.
- **Summary:** Viewed as an inverse problem, each stage of the pipeline has more freedom than the data it is fit to can constrain, and several "locks" are one-number calibrations of a curve whose shape is the actual content. The paper itself shows the symptom: six density laws give the same gradient (−0.00106 to −0.00109 Gpc⁻¹) because only ρ′/ρ at D is constrained, and the cluster data are matched equally well by laws with very different global behaviour.

## Strengths
- **S1** The level-set lock is mathematically clean: one constraint fixes one scale per law; all six reproduce (M11–M13).
- **S2** The axisymmetric curve A_F is genuinely the CoE forward map, not a disguised fit (M16) — a fair test object.
- **S3** Good use of a Δχ² profile (Fig. 9) to show the L–D degeneracy.

## Weaknesses
- **W1 [MAJOR] Identifiability.** Location: §4.1, Figs. 8–9. The lock fixes only the local log-gradient; the cluster anisotropy fixes mainly a dipole + quadrupole phase. A 2-parameter family (D, L) is then being constrained by ~2–4 effective numbers. Fix: report the Fisher matrix / posterior for (D, L, shape, l₀) using a proper covariance and show which combinations are constrained.
- **W2 [MAJOR] Unspecified forward operator.** Location: Eq. 3 (f, p(s) unspecified, self-flagged). The core map from density to A cannot be reproduced from the paper (M9). Fix: give f and p(s) explicitly with code.
- **W3 [MAJOR] Free-form source in the CMB stage.** Location: Eq. 19, Table 2. A 12-parameter P_δ with two oscillatory terms and a free projection centre is a flexible basis; reaching r²=0.976 shows the basis spans the TT shape, not that the model predicts it. Fix: a prediction test — fit ℓ<800, predict 800–2500 (and TE/EE) blind.
- **W4 [MAJOR] Two sky orientations.** Location: §4.6 vs §4.7.2. Fitting uses +80.8° roll; anomaly comparison uses 15° lift + 45° rotation, unreconciled. Fix: one orientation, fixed before the comparison.
- **W5 [MINOR] D=42 vs 58.2** in the web modulation (M17) — the web is built for a different universe than the one the cluster fit selects.
- **W6 [MINOR] Pass/fail rule mislabelled** (R=8 PASS with lowl_over 2.99 > 2.5; R=14 FAIL unstated) (M25).

## Alternative perspective
The same pipeline run with a pure ΛCDM + bulk-flow (or kinematic dipole + quadrupole) template on the cluster data would likely do at least as well (cf. M16: 3-harmonic χ²=3.17). The paper would be stronger as a demonstration that "an off-centre matter geometry is *not excluded* by cluster scaling data at the level of a low-order template", with the CMB stage separated out as an open problem.

## Questions
1. How many effective degrees of freedom does the A(n̂) curve have after the 90° cone smoothing?
2. Can the CoE model predict anything not used in calibration (e.g. the H₀ anisotropy amplitude in a different cluster sample)?
