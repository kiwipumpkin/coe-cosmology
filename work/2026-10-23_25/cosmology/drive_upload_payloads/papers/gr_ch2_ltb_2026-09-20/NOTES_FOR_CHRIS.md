# Notes for Chris — GR Ch.2 LTB shot (2026-09-20)

Plain speech.

## What worked

- Got a real LaTeX chapter on disk with LTB metric, decreasing-\(\rho\) toy \(M(r)\), off-center geometry, sheet-AF pipeline, Migkas comparison, and a CMB kill-shot section.
- Used the digitised Migkas CSV; did not invent \(A(\ell)\) points.
- Two required figures exist and are in the PDF:
  - toy \(A(\ell)\) vs Migkas
  - CMB dipole vs offset \(d\)
- Citations: Alnes, Kodama/Saito, Bertschinger, Smoller–Temple, Migkas, Heinesen, Saadeh (Bianchi = what we are *not* doing). Acknowledgments name Rindler and Ackerman.

## What is still Newtonian scaffolding

- The “LTB” numerics are **not** full null geodesics through an LTB metric. They evaluate an FLRW-shell \(H(r)\) built from enclosed mass \(M(r)\) on an observer-centered source circle, then map \(A\propto H_\mathrm{obs}/H\) exactly like Chapter 1 sheet AF.
- So the angular template is still mostly the Newtonian enclosed-mass idea wearing GR clothes. A finished chapter needs past-directed null rays, \(d_L(z,\hat n)\), and a real temperature map.
- Amp×sin (and 2ℓ if used later) are still diagnostic knobs, not Einstein-equation outputs. In this toy Amp×sin barely helps.

## Numbers (from `toy_results.json`)

| Item | Value |
|------|-------|
| Geometry + phase affine \(R^2\) vs Migkas \(A\) | **≈ 0.430** |
| + Amp×sin diagnostic \(R^2\) | **≈ 0.494** |
| Best phase | \(84^\circ\) |
| \(D/L\), \(s/D\) | 6, 0.5 |
| CMB: \(d_\max\) for observed dipole (this profile) | **≈ 2.5 Mpc** |
| Dipole at \(d=0.5D\sim1500\) Mpc (linear toy) | **≈ 0.73** (fatal) |
| Alnes-scaled at same \(d\) | **≈ 0.14** (also fatal) |
| Alnes/Kodama void ceiling class | **~15 Mpc** |

Do not compare the 0.43 \(R^2\) to Chapter 1’s 0.54–0.88 ladder as a win/loss of GR vs Newton — different proxy, weaker angular content. Report it straight.

## CMB verdict from the toy

**Spherical beyond-horizon CoE + off-center us at Newtonian \(d/D\sim0.5\) loses to the CMB dipole.** Horizon-scale offsets produce dipoles orders of magnitude too big. Surviving would need tiny offset (\(\mathcal{O}(10)\) Mpc), broken spherical symmetry with care, or fine-tuned cancellation — none of which this shot demonstrates.

## Open next

1. Real LTB geodesic \(d_L\) and \(\delta T/T\).
2. Bertschinger self-similar blast \(M(r,t)\) instead of static exponential.
3. Smoller–Temple matching for beyond-horizon ontology — after the local observable pipeline is honest.
4. Do **not** pivot to Bianchi as a loophole; Saadeh already crushed that class.
