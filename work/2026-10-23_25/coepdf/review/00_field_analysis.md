# Phase 0 — Field Analysis & Reviewer Configuration (field_analyst_agent)

## Paper Basic Information
- **Title:** Off-centre CoE cosmology: the mathematical workflow from geometry to the Planck TT spectrum (rev. 5, compiled 24 Oct 2026)
- **Source:** /workspace/coepdf/main_rev5.tex (431 lines), PDF /workspace/coe_math_workflow_rev5.pdf
- **Type:** methods/working-notes manuscript; phenomenological model + fits to three data sets (Migkas L_X–T anisotropy, CatWISE dipole, Planck 2018 TT).
- **Maturity:** working draft. Self-flags (`[unverified: …]`) present; references incomplete (no years/journals for most entries); two sky orientations unreconciled.

## Field Analysis
| Dimension | Assessment |
|---|---|
| Primary discipline | Physical cosmology (inhomogeneous / off-centre cosmologies, CMB phenomenology) |
| Secondary | Observational X-ray cluster cosmology; large-scale-structure statistics; inverse problems / statistical model checking |
| Paradigm | Quantitative, model-fitting, exploratory |
| Methodology type | Phenomenological forward modelling (Newtonian + FRW thermal history), χ² minimisation, harmonic residual analysis, Monte-Carlo random webs |
| Target tier | As written: arXiv preprint / working notes. With the issues below fixed: MNRAS, PRD, JCAP (all routinely publish off-centre LTB and anisotropy tests) |

## Recommended Target Journals (Top 3)
1. Phys. Rev. D — standard home for LTB/off-centre-observer CMB constraints (Alnes & Amarzguioui 2006; Foreman et al. 2010).
2. JCAP — anisotropy/dipole-tension papers.
3. MNRAS — Migkas-type cluster anisotropy follow-ups.

## Reviewer Configuration Cards
### Card #1 — Journal-Fit Reviewer (EIC)
Editor of a PRD/JCAP-type journal handling non-standard cosmologies. Focus: is the claim novel against the LTB-void literature, is it framed honestly, is it a paper or notes.
### Card #2 — Reviewer 1 (Methodology)
Cosmological statistician (likelihoods, correlated binned data, model comparison, reproducibility). Focus: χ² on overlapping sliding cones, dof accounting, look-elsewhere in the 11/12 anomaly test, ensemble-mean vs realisation, code-to-number traceability.
### Card #3 — Reviewer 2 (Domain)
CMB / inhomogeneous-cosmology theorist (LTB, Sachs–Wolfe, lensing, kinetic SZ). Focus: CMB isotropy bounds on off-centre observers, GR validity at D ≈ 58 Gpc, consistency of using ΛCDM inputs (Ω_m incl. dark matter, Eisenstein–Hu CDM transfer) inside a "no dark matter" model, literature coverage.
### Card #4 — Reviewer 3 (Perspective)
Applied mathematician in inverse problems / geophysical tomography. Focus: identifiability of m(r), degeneracy, regularisation, whether residual-harmonic peeling is diagnostic or just flexible curve fitting.
### Seat #5 — Devil's Advocate (fixed)
No card. Attacks core thesis and logic chain.

## Review Strategy Recommendations
- Methodology and DA seats carry the most weight (statistics on Migkas bins and the CMB stage are where the claims rest).
- A separate mathematics check (re-derivation + numerical spot check against box code/data) is attached as Appendix M.
- Provenance note: dispatch of independent subagents was unavailable; the five Phase-1 reports were written sequentially by one executor, each written to file before the next began. Role separation is therefore procedural, not independent error processes (as the skill itself states).
- Calibration status: NOT_CALIBRATED.
