# Module 3 — Distance measures and redshift

**Audience.** Chris: chemistry → applied physics.

**Stance.** No cheerleading of ΛCDM (Lambda cold dark matter) or of the Center-of-Expansion (CoE) sketch.

### Acronyms

| Acronym | Meaning |
|---------|---------|
| **CP** | Cosmological principle |
| **CoE** | Center of Expansion |
| **FLRW** | Friedmann–Lemaître–Robertson–Walker |
| **CMB** | Cosmic microwave background |
| **SN** / **SNIa** | Supernova / Type Ia supernova |
| **CM** | Center of mass |

**Spine.** \(D\) = CoE→us; \(r\) = vision shell; \(R_{\mathrm{CMB}}\) = last-scattering / **CMB** shell; decreasing \(\rho\); Gauss/shell.

---

## 1. Goal

Master Hogg’s distance-ladder language (comoving, luminosity, angular-diameter) and see how an **off-center** observer on a CoE-centered vision shell confuses directional distance and inferred \(H\).

---

## 2. Chem / applied-physics bridge

**Lab frame vs center-of-mass (CM) frame:** the same collision looks different once you boost. Cosmological redshift mixes expansion, peculiar velocity, and path-dependent gravity. Frame choice matters — especially when you are not at the symmetry center (\(D\neq 0\)).

---

## 3. Distance measures (FLRW language first)

In a **Friedmann–Lemaître–Robertson–Walker (FLRW)** cosmology (from the **cosmological principle (CP)**), standard distances include:

- **Comoving distance** — “ruler on the expanding grid.”
- **Luminosity distance** \(d_L\) — from flux and intrinsic luminosity (\(d_L\) enters supernova (**SN**) Hubble diagrams).
- **Angular-diameter distance** \(d_A\) — from physical size and angular size.

In FLRW with standard assumptions they are related by redshift factors (Etherington’s distance-duality in transparent GR). Once expansion history is nontrivial, \(d_L\) and \(d_A\) disagree numerically even though duality holds — that is history, not a bug.

**Hogg (1999)** is the clean reference: David W. Hogg, “Distance measures in cosmology,” arXiv:astro-ph/9905116.

---

## 4. Project thread — off-center shell → directional distance / \(H\)

Sit at distance \(D\) from the CoE. Sources on a shell of radius \(r\) from the CoE are **not** all at the same distance from you:

- Toward the CoE: path length roughly related to \(|r - D|\) (geometry depends on whether \(r > D\)).
- Away from the CoE: roughly \(r + D\).
- Off-axis: intermediate.

If you naively convert redshift → distance with an FLRW formula that assumes isotropy about *you*, a real variation of path length / expansion along different lines of sight can appear as a **sky dipole in inferred \(H\)** or in \(d_L\). That is the Module 3 spine item linking to Migkas-style directional claims (Module 5) and to the CMB path asymmetry \(R_{\mathrm{CMB}}\pm D\) (Module 6).

Figure (shared machine): `/workspace/cosmology/figures/coe_geometry_cmb_dipole_schematic.png` — the same \(R\pm D\) intuition that budgets the CMB dipole also budgets directional path length to a vision shell.

**Honesty.** Until you integrate past-directed null geodesics in a real metric, “\(H\) on the sky” from enclosed-mass Newtonian shells is a **diagnostic template**, not a finished GR observable (see Module 1 / NOTES_FOR_CHRIS).

---

## 5. Readings

- David W. Hogg, arXiv:astro-ph/9905116 (1999; rev. 2000).
- Barbara Ryden, *Introduction to Cosmology*, distance and redshift chapters.
- Steven Weinberg, *Cosmology* (2008), luminosity distance.
- **Project thread:** `/workspace/cosmology/figures/coe_geometry_cmb_dipole_schematic.png`.

---

## 6. Exercises

1. In words: why can luminosity distance and angular-diameter distance disagree once expansion history is nontrivial?

2. **Project thread: Off-center shell → directional \(H\).** You sit at distance \(D\) from the CoE. Sources on a shell of radius \(r\) from the CoE are at different distances from *you*. Which observable can look like a sky dipole in \(H\)? Sketch with \(r\) and \(D\) labeled.

3. Name one quantity that is frame-dependent (lab-frame velocity) and one closer to an invariant (proper interval along a null geodesic). Map “lab offset from blast center” to \(D\).

---

## 7. Takeaway

Distances are defined carefully in FLRW; off-center geometry makes “the” distance to redshift \(z\) direction-dependent. Module 4 separates bulk flow from true directional expansion.
