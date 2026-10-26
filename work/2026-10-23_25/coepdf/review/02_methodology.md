# Peer Review Report — Reviewer 1 (Methodology)

## Reviewer Information
- **Role:** Methodology (card #2). **Identity:** cosmological statistician. **Focus:** likelihoods, correlated data, dof, model comparison, reproducibility.

## Overall Assessment
- **Recommendation:** Major Revision. **Confidence:** 4/5. **Calibration:** NOT_CALIBRATED.
- **Summary:** The arithmetic is largely right (see Appendix M), but the inferential layer is not. The cluster χ² is computed on 71 heavily overlapping sliding-cone bins as if independent; the 11/12 anomaly score has no null distribution; the CMB stage compares an ensemble-mean spectrum to one sky while labelling it a realisation; and the method's first step (run9 shear) does not reproduce.

## Strengths
- **S1** Parameters and χ² values in Table 2 reproduce exactly from `acoustic_half_fixed.json` (M36).
- **S2** Cosmic-variance treatment (Eq. 8, χ²_CV) is correct and its limits stated.
- **S3** Limber validity is quantified against the exact projection (acoustic_limber_check*.json).

## Weaknesses
- **W1 [CRITICAL] Correlated Migkas bins treated as independent.** Location: §3.2, §4.1 "Fit", Abstract. Evidence: Migkas Fig. 6 cones have Δl=90° and are sampled every 5°, so each cluster contributes to ~18 adjacent bins; effective independent bins ≈ 360/90 = 4. Under the independence assumption the null χ²=45.14/71 has P=0.993 (the null is "too good"), and fits give χ²/dof ≈ 0.11–0.17 (M15). Consequence: χ² differences 45→7.67 have no probabilistic meaning; with ~4 effective points, a 2–4 parameter model fitting them well is expected. Fix: refit the 313-cluster catalogue (migkas2020_tableC1.csv is on the box) directly: per-cluster L_X–T likelihood with a direction-dependent A(n̂), or build the full 71×71 bin covariance by bootstrap over clusters and use χ²=rᵀC⁻¹r.
- **W2 [MAJOR] No comparison against generic low-order models.** Location: §4.1, Fig. 7. Evidence: a free 3-harmonic Fourier series (7 params) gives χ²=3.17 on the same 71 bins; the paper's own best "variant" (7.67) is said to equal a "reference sinusoid". Fix: report AIC/BIC or cross-validated likelihood of every CoE law against n=1,2,3 Fourier baselines; the CoE geometry must beat equally flexible phenomenology to count as evidence.
- **W3 [MAJOR] 11/12 sky-anomaly score lacks a null and has researcher degrees of freedom.** Location: §4.7.2, Fig. 22, Abstract. Evidence: orientation (15° lift, 45° rotation), smoothing (±25°), harmonic order (n=1–3), exclusion regions and A/B assignment are all choices; antipodal spots (dipole poles, bulk-flow toward/away) are anti-correlated by construction; one spot sits on zero (+0.02 %). Naïve binomial P=0.0032 is not valid. Fix: Monte-Carlo the score over random rotations of the residual map (and random spot assignments), report the fraction scoring ≥11/12; pre-register the orientation.
- **W4 [MAJOR] Ensemble mean vs realisation.** Location: Fig. 18 caption, Table 2. Evidence: χ²=3,787 is the ensemble-mean spectrum; the seed-42 one-sky realisation gives 4,601 (2–2500) and 148 (2–30) (M37). Fix: label correctly; compare the observed sky to the distribution of realisations.
- **W5 [MAJOR] First fit step not reproduced and not reported.** Location: absent (method §4.4/4.5 never presents the total-shear run9 step). Evidence: rerun gives 272,638 bodies and pattern/noise 48.9 vs recorded 272,175 / 149.9 (M26). Fix: pin code versions/seeds, archive inputs, report the shear step and its sensitivity to seed and thinning.
- **W6 [MAJOR] Parameter counting and boundary solutions.** 12 CMB parameters + D, l₀, shape, tilt, gradient, roll, rotation; "half" fit has s_w on its bound (0.050) (M39). The "fixed" sky placement (tilt/gradient/roll) is described as adding "no morphological degree of freedom", yet the 80.8° roll and 45° rotation are chosen values. Fix: count them as parameters; report bounds hit.
- **W7 [MINOR] Error convention.** 45.14 uses the mean of asymmetric errors; Eq. 9's side rule gives 43.4 (M14). State it.
- **W8 [MINOR] Data file.** §3.3 says binned file for ℓ≥30, fit uses 2,499 unbinned points (M36).
- **W9 [MINOR] Many results "taken from earlier figure set; script not re-run"** (harmonic amplitudes, ℓ=50 step, Δℓ≈1384). Re-run or remove.

## Questions for Authors
1. What is the bin-to-bin correlation matrix of the digitised Fig. 6?
2. What fraction of random map rotations give ≥11/12?
3. Which code revision produced the repo's 149.9?

## Criterion-Bound Judgements
| Criterion | Judgement |
|---|---|
| Design validity | Weak (inference on correlated bins) |
| Statistical reporting | Inadequate (no significance, no covariance, no model-selection penalty) |
| Reproducibility | Partial (CMB numbers reproduce from JSON; shear step does not; several scripts missing) |
