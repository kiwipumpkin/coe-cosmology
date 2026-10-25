# Synthesis note — ST + CoE paper (2026-09-19)

**For:** Chris  
**Output:** `papers/st_coe_synthesis_2026-09-19/main.pdf`  
**Author line:** C. A. (Chris), matching `new_draft.tex` voice.

## Intellectual arc (how this was synthesized)

1. **Ontology first.** The unfinished draft and Phase A derivation already had the binding claim: a *real* Center of Expansion beyond the horizon, not a ΛCDM dipole patch. The synthesis keeps that as the thesis, not a mild tweak.

2. **Kinematic honesty.** Phase A’s homology theorem (`v = H r` ⇒ isotropic `H_obs`) is stated up front. Anisotropy is attributed only to Gauss gravity with *decreasing* ρ — uniform-ball `g ∝ −r` is explicitly rejected. Numbers for κ ~ O(0.1–1) vs M_crit come from the absurdity / enclosed-mass notes, not invention.

3. **Bolt-on → ideal integration.** `integrated_gravity.py` supplies the closed α_eff results (frozen → 2, drag → 2/3). The paper’s structural punchline is Chris’s: integration rescales amplitude only; expand vs collapse flips A(φ). Sheet AF uses `A ∝ H_obs/H`; invsq is mentioned as the trough/peak inverter.

4. **Empirical ladder.** Fit values are taken from the same chain as `best_fit_2ell_multiplier.npz` and the ideal_rotation / Amp–φ CSVs:
   - geometry alone R² ≈ 0.543  
   - Amp×sin (Amp=−0.06, φ≈0.20 rad) R² ≈ 0.746  
   - + 2ℓ (ε≈−0.043, ψ≈176°) R² ≈ 0.875, RMSE 0.038→0.029  
   - top geometry M/M0=3, d/D=0.5, L≈0.169  

5. **Knob honesty + CMB adversary.** Amp×sin and 2ℓ are labeled diagnostic; nesting dolls speculative. Related-work section compresses `PAPER_DIGEST.md` (Bertschinger, Smoller–Temple, Alnes/Kodama, Saadeh ≠ inhomogeneous CoE, Migkas Fig. 6). Acknowledgments restore Rindler & Ackerman as genesis.

## Sources used (not re-derived)

- `/workspace/cosmology/phase_a_derivation.md`, `new_draft.tex`, `procedure_methods.tex`
- `density_local_match/run_density_local_match.py` (`H_obs_circle`, `af_an_from_physics`)
- `integrated_gravity.py`, `absurdity_check.md`
- `map_compare/amp_mean_param_noise/*` figures + NPZ/CSVs
- `literature/related/PAPER_DIGEST.md`

## What this paper is not

A discovery claim, a GR metric, or a CMB-compatible construction. It is a peer-review-toned synthesis of the Newtonian program with explicit R² ladder and falsifiers.
