# Our sphere over/underdense, with giant partners outside it

Chris: "If our universe is on avg over dense or under dense, it'll have associated under or over dense partners. Try those type of extreme masses or voids next."

## Setup
- **Our sphere:** the 14.3 Gpc visible sphere gets a uniform density contrast δ_us, fitted between −0.3 and +0.3. Mass per unit δ_us is 1.04e54 kg.
  - A perfectly even sphere around us has no effect on the ring. Only the slight density slope toward the CoE leaves a pattern: 515× floor per unit δ_us, almost all of it a dipole (m = 1). It is computed numerically in `sphere_pattern.py` and converges to under 1%.
- **Partners:** 1–3 spheres with radius 1–15 Gpc, with centres up to 60 Gpc from us. Each must lie wholly outside our sphere, outside the CoE interior and inside the shell's outer edge.
  - "Overdense us": partners are voids (δ −0.9 to 0).
  - "Underdense us": partners are overdensities (δ 0 to +2).
  - "Free": either sign.
- **Net mass is exactly zero:** sphere excess + partners + any ordinary objects must sum to zero (|net| ≤ 1e51 kg). There is no smooth background compensation, so the partners must do the compensating.
- **Ordinary objects (optional):** the usual masses (≤ 2e52 kg) and deep voids (≤ 4 Gpc), as in `../ss_core.py`.
- **Starts:** 24 per partner-only setup and 16 or more per combined setup, mixing warm starts from the old 4-mass and 4+4 fits with cold starts.

## Results (`a_best_by_setup_logy.png`)
- **Partners alone don't help.** The best result is 106× (3 free partners), against 237× at the start. The fit keeps δ_us ≈ 0 and shrinks the "partners" to ordinary 2–4 Gpc objects about 17–22 Gpc away. These do worse than the earlier 2–3 ordinary objects (92–121×).
  - The only giant object in this group is in "overdense, 3 partners" (153×). There δ_us sits at the +0.3 limit and is balanced by a 9.8 Gpc void at ℓ 195°, b −61°.
- **With 4 masses + 4 voids, the best is 1.23×** (free and underdense setups find the same layout).
  - δ_us = −0.025, i.e. our sphere at almost exactly average density.
  - The 2 partners are compact overdensities: R 3.2 Gpc at ℓ 297°, b +23°, 28 Gpc out (2.3e52 kg), and R 2.4 Gpc at ℓ 159°, b −24°, 19 Gpc out (1.0e52 kg). In effect these are 2 more masses in the usual bands.
  - That gives 10 objects at 1.23×, the same as the old 10-object result (1.26×). The partners don't reduce the number of objects needed.
- **Best layout that really uses a giant partner** (`d_best_giant_partner_map.png`): overdense us, 2 partners + 4 masses + 4 voids, **1.57×**.
  - δ_us = +0.21, a sphere excess of 2.2e53 kg.
  - The main partner is a giant void: R 9.2 Gpc at δ −0.9 (the depth limit), ℓ 190° b +66°, centre 25.7 Gpc from us, touching the CoE interior edge. Its deficit of −2.5e53 kg outweighs the sphere excess plus the ordinary masses.
  - It does slightly worse than the near-average-density layout.
- **Mass bands:** the ordinary masses again sit in the usual bands (ℓ ≈ 75–90°, 160°, 220–230°, 300–315°).

## Answers
- **Preferred δ_us:** close to 0 (−0.02 in the best fits).
  - Underdense us gets to 1.23× but uses only δ ≈ −0.02.
  - Overdense us with a giant void gets to 1.57× at δ +0.21.
  - Large |δ_us| is mainly something the fit has to cancel: its dipole is 730× per unit δ, while the target dipole is only 18×.
- **Partner sizes:** the data never asks for 10–15 Gpc partners unless δ_us is pushed large. The radius never reached 15 Gpc. Contrasts often sit at the −0.9 / +2 limits, and the 4 ordinary masses often sit at the 2e52 cap.
- **Mass compensation:** holds by construction (|net| ≤ 1e51 kg in every fit). In the overdense setups, the giant void has to absorb the sphere excess and the masses' 8e52 kg.

## Cut for time
- A wider δ_us range (up to +0.6) was prepared but not run.
- Some "free" combined starts were duplicated by a parallel helper run (`raw_b/`). These were merged, not lost.

Files:
- Code: `ep_core.py`, `sphere_pattern.py`, `run_ep.py`, `run_ep_b.py`, `report_ep.py`
- Results: `raw/`, `raw_b/`, `ep_results.json`, `report_ep_out.txt`
- Charts: `d_best_partner_map.png` (overall best), `d_best_giant_partner_map.png`
