# Editorial Decision

**Manuscript:** Off-centre CoE cosmology: mathematical workflow (main_rev5.tex)
**Decision:** **MAJOR REVISION** (as a journal article it would be Reject-and-resubmit; the three CRITICAL items must be resolved before re-review)
**Calibration:** NOT_CALIBRATED. Provenance: all five passes were written sequentially by one executor following the ARS academic-paper-reviewer workflow; they are not independent human reviews.

## Blocking Issues (3)
1. CMB isotropy for a 58-Gpc off-centre observer not demonstrated (DA-2, R2-W1). [CONSENSUS-3: DA, R2, EIC]
2. Cluster χ² on correlated bins; null not rejected; generic harmonics fit better (DA-1, R1-W1/W2, R3-W1). [CONSENSUS-4]
3. Inconsistent cosmological background (H₀≈47 vs 70 vs ~145; ΛCDM dark-matter inputs) (DA-3, R2-W2/W3). [CONSENSUS-3]

## Reviewer Summary
| Reviewer | Recommendation | Confidence |
|---|---|---|
| EIC (journal fit) | Major Revision | 4 |
| R1 Methodology | Major Revision | 4 |
| R2 Domain | Major Revision (Reject-leaning for CMB) | 4 |
| R3 Perspective | Major Revision | 3 |
| Devil's Advocate | Reject in current form | 4 |

[SPLIT] Severity of the CMB acoustic fit: EIC/R3 treat it as a legitimate "can it mimic" exercise if reframed; R2/DA treat it as evidence against the model. Editor: reframe, and state Δχ²≈1,229 as disfavouring.

## Top 10 issues (ranked)
| # | Severity | Location | Problem | Concrete fix |
|---|---|---|---|---|
| 1 | Critical | §4.3, §4.5 | Residual geometric CMB multipole after 95.9 % cancellation ≈ 4e-4–2e-3 ≫ 1e-5; off-centre LTB bounds (observer within ~tens of Mpc of centre) not addressed | Compute ℓ=1–3 ΔT/T of the full final sky; confront Alnes & Amarzguioui 2006, Foreman et al. 2010 |
| 2 | Critical | §3.2, §4.1, Abstract | 71 overlapping 90° cones at 5° steps (~4 independent) fit as independent; χ²/dof ≈ 0.11–0.17; null P=0.993 | Per-cluster likelihood on the 313-cluster catalogue or bootstrap 71×71 covariance |
| 3 | Critical | Eq. 2, Eq. 5, §4.4, Table 1 | t₀=13.8 Gyr w/o Λ ⇒ H≈47; lock uses H₀=70; R=4 Gpc ⇒ H₀≈145; web uses Ω_m=0.31 + CDM transfer | One baryon+radiation background; derive a_eq, R, growth, σ₈ from it |
| 4 | Major | §4.1, Fig. 7 | 3-harmonic Fourier (χ²=3.17) beats all CoE laws (7.67–11.6) | Model comparison (AIC/BIC/CV) vs n=1–3 templates and a bulk-flow model |
| 5 | Major | §4.7.2, Abstract | 11/12 sign score: post-hoc orientation/smoothing/assignment, correlated antipodal spots, no null (naïve binomial 0.0032 invalid) | Rotation/assignment Monte Carlo; one pre-fixed orientation |
| 6 | Major | Eq. 19, Table 2 | Acoustic peaks inserted via cos(k r_s), r_s free; χ² 3,787 vs 2,558 (χ²/dof 1.52 vs 1.03) | Reframe as mimicry test; blind high-ℓ / TE/EE prediction |
| 7 | Major | §4.5; run9 | Tidal shear/lensing omitted from TT stage; shear step not reported and not reproduced (48.9 vs 149.9) | Add lensing; pin seeds/code; report and reproduce shear step |
| 8 | Major | Eq. 4 | CatWISE dipole (kinematic+clustering of z≈1–2 sources) equated to density contrast at D±14.3 Gpc | Integrate model density+velocity over CatWISE dN/dz |
| 9 | Major | Fig. 18 caption, Table 2 | 3,787 is the ensemble mean; seed-42 realisation gives 4,601 (2–30: 148) | Relabel; compare data with realisation distribution |
| 10 | Major | §1, refs; Abstract | Motivating SN acceleration never computed; LTB/kSZ/μ-distortion literature absent | Compute D_L(z,n̂) vs Pantheon+; add related work |

Further minor: Eq. 19 example 2.2e-9 → 2.06e-9; Planck peak lists disagree (224/535/821/1122/1412 vs 228/544/823/1121/1417); D=42 vs 58.2 in web; R=8 labelled PASS despite lowl_over 2.99>2.5 (and R=14 FAIL unstated); s_w on bound; ℓ notation clash; χ² error convention; unbinned vs binned data; δ rounding (L 920.3 vs 922.65); incomplete references; f, p(s), K_G unspecified.

## DA-CRITICAL disposition
DA-1, DA-2, DA-3 accepted as CRITICAL (data–conclusion mismatch; foundation collapse; logic-chain break). DA-4–6 accepted as Major.

## Rationale
The arithmetic is sound (≈35 of 45 checks verify; see Appendix M) and the manuscript is unusually candid. But the inferential claims (cluster improvement, sky-anomaly agreement) are not supported by valid statistics, the model is not internally consistent, and the most stringent test (CMB isotropy) is unaddressed. These are fixable in principle but require new analysis, not edits.

## Roadmap
1. Fix the background (Issue 3). 2. Full-sky low-ℓ CMB (Issue 1). 3. Cluster likelihood + model comparison (2, 4). 4. Anomaly null (5). 5. Reproduce shear, add lensing (7). 6. Reframe TT fit and relabel (6, 9). 7. Minor list.
