# Peer Review Report — Reviewer 2 (Domain)

## Reviewer Information
- **Role:** Domain (card #3). **Identity:** CMB / inhomogeneous-cosmology theorist. **Focus:** physical consistency, CMB isotropy, literature.

## Overall Assessment
- **Recommendation:** Major Revision (bordering Reject for the CMB section as a physical model). **Confidence:** 4/5. **Calibration:** NOT_CALIBRATED.
- **Summary:** The model mixes three incompatible frameworks: Newtonian kinematics at distances where recession is ≫c, an FRW thermal history normalised to an age that implies H₀≈47, and ΛCDM inputs (Ω_m=0.31 including dark matter, σ₈, Eisenstein–Hu CDM transfer, a_eq=1/3400) in a model stated to have no dark matter. The CMB spectrum is then reproduced by a phenomenological P_δ with an acoustic-like oscillation and a second damped mode inserted by hand, so the acoustic peaks are put in, not predicted.

## Strengths
- **S1** Correct statement that geometry alone gives a percent-level CMB multipole (§4.3) — consistent with LTB literature.
- **S2** Correct linear SW/ISW/Limber machinery (M28–M31).
- **S3** The paper recognises that the moving-lens regime is superluminal and treats it as an assumption.

## Weaknesses
- **W1 [CRITICAL] CMB isotropy for an observer 58 Gpc off-centre.** Location: §4.3. Evidence: even after the stated 95.9 % cancellation, a 1–5 % geometric multipole leaves 4e-4–2e-3 (M43), 40–200× the observed ~1e-5. Off-centre LTB studies (Alnes & Amarzguioui 2006; Foreman et al. 2010; Quartin & Amendola 2010) restrict the observer to within tens of Mpc of the centre of a Gpc structure using dipole/quadrupole/octopole. Fix: compute the full dipole–quadrupole–octopole of the final model's CMB sky (geometry + Doppler + gravitational redshift + web), not just the stochastic part's power spectrum, and confront the off-centre bounds explicitly.
- **W2 [CRITICAL] Internal inconsistency of cosmological inputs.** Location: §2.2 Eq. 2, §4.1 CatWISE lock, §4.4 ingredients. Evidence: t₀=13.8 Gyr without Λ requires H₀≈47 km/s/Mpc; the lock uses H₀=70, Ω_m=0.3; the web uses Ω_m=0.31, Ω_b=0.049 (i.e. 84 % non-baryonic matter) and a CDM transfer function (M5, M22). Fix: pick one background with only baryons + radiation (Ω_b-only), derive a_eq, t(a), growth, transfer function and σ₈ from it, and redo the web.
- **W3 [MAJOR] Recombination radius not derived.** Location: Table 1, §2.1, §4.5. Evidence: R=4 Gpc requires H₀≈145 in a matter-only model; the model's own value would be 8–12 Gpc (M6). Fix: derive R from the CoE kinematics/D(t) and keep it fixed in the TT fit.
- **W4 [MAJOR] Acoustic peaks inserted by parametrisation.** Location: Eq. 19. Evidence: A_w cos(k r_s+φ) and B e^{-k/k_B}cos(k r₂+φ₂) supply the peak structure; r_s=79.5 Mpc is free. In ΛCDM the peaks, their odd/even heights (baryon loading) and Silk damping are predicted from Ω_b h², Ω_c h². Fix: state clearly that the TT fit tests whether a parametrised source can mimic the spectrum, not whether the CoE model predicts it; or derive P_δ from baryon-photon physics in the CoE background.
- **W5 [MAJOR] Tidal shear/lensing absent from TT stage.** Location: §4.5. Evidence: fwd.py uses SW+ISW+moving lens only (M44); the shear step (run9) is separate and unreported. Planck measures CMB lensing (C_ℓ^φφ) and peak smoothing consistent with ΛCDM at the few-% level. Fix: add lensing of the primary spectrum and predict C_ℓ^φφ; compare with Planck lensing.
- **W6 [MAJOR] Missing observables.** Polarisation (TE/EE), kSZ (Zibin & Moss 2011 constraints on large-scale radial flows, especially with a ≥4c moving-lens flow), μ/y spectral distortions (Caldwell & Stebbins 2008), BAO, SN D_L(z). Fix: at least discuss kSZ and spectral distortions, which are the decisive LTB tests.
- **W7 [MAJOR] CatWISE lock physics.** Location: Eq. 4. Evidence: the number-count dipole is a kinematic + clustering signal of sources at z≈1–2, integrated along the line of sight, not ρ(D−d)−ρ(D+d) at 14 Gpc (M10). Fix: model the predicted CatWISE dipole by integrating the CoE density and velocity over the quasar dN/dz with α, x from Secrest et al.
- **W8 [MINOR]** Follow-up dipole work (Secrest et al. 2022, NVSS+CatWISE) and Migkas et al. 2021 (bulk-flow interpretation of the same anisotropy) are relevant and absent. [citations from memory — verify details]

## Questions for Authors
1. What is ΔT/T at ℓ=1,2,3 for the full final model sky?
2. What fixes r_s=79.5 Mpc physically?
3. How is a superluminal Newtonian moving-lens kernel justified versus a relativistic LTB ray-trace?

## Criterion-Bound Judgements
| Criterion | Judgement |
|---|---|
| Theoretical consistency | Fails (mixed frameworks) |
| Literature coverage | Inadequate |
| Contribution | Exploratory; no predictive CMB result |
