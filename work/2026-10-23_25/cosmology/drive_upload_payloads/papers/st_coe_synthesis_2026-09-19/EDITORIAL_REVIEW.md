# Editorial review — ST + CoE synthesis (2026-09-19)

**Paper:** `main.tex` → `main.pdf` (14 pages)  
**Editor pass:** submission-readiness audit (figures, math, claims, prose, refs)

---

## Critical fixes made

1. **Expand/collapse figure (Chris-flagged).** Removed `expand_vs_collapse_vs_migkas.png` from the paper. That file’s bottom panel is Migkas **B** (`B/B_all`), wrongly implying a Migkas-B comparison. Replaced with `figures/A_vs_angle_ideal_collapse.png` (copied from `map_compare/amp_mean_param_noise/`). New caption: **top** = collapse vs expand mean-renorm \(A(\varphi)\) under ideal drag (Migkas for scale); **bottom** = sign flip of \(f = H_{\mathrm{obs}}/H_{\mathrm{bg}} - 1\) — **not B, not Migkas B**.

2. **`migkas_bars_rotation_overlay.png` honesty.** Bottom panel *is* truly Migkas \(B\). Caption updated: A top / B bottom; B marked ancillary and **not** used in the \(R^2\) ladder; note that panel phase \(\phi\approx0.35\) sits on the Amp\(=-0.06\) ridge, while the ladder quote is \(\phi\approx0.20\).

3. **Caption / plot mismatches elsewhere.**
   - `A_vs_angle_ideal_integration_sheet_AF.png`: caption now states AF \(= A_{\mathrm{hratio}}\propto H_{\mathrm{obs}}/H\) and the illustrative \(\kappa\approx0.231\), \(d/D=0.1\) slice.
   - `A_vs_angle_ideal_plus_chris_rotation.png`: caption clarifies sheet-style AF + \(\phi\approx0.35\) (not the Table ladder peak \(R^2\approx0.746\) at \(\phi\approx0.20\)).
   - `dipole_vs_quadrupole_examples.png`: caption clarifies schematic multipole coeffs \(A\)/\(B\) ≠ Migkas strip \(A\)/\(B\).
   - Geometry-grid caption now quotes \(M/M_0=3\), \(d/D=0.5\), \(L\approx0.169\).

4. **Claims / prose discipline.** Abstract and conclusions tightened; informal “smoked the CP” → “falsified the cosmological principle”; Amp\(\times\sin\) / \(2\ell\) remain diagnostic; nesting dolls speculative; CMB adversary (Alnes/Kodama \(\sim15\,\mathrm{Mpc}\)) retained. Author line `C.\ A.`; Acknowledgments keep Wolfgang Rindler and Luke Ackerman.

5. **Build.** `pdflatex` ×2 succeeded; unused `expand_vs_collapse_vs_migkas.png` left on disk but not included. Hyperref math-in-bookmark warnings fixed via `\texorpdfstring`.

---

## Consistency checklist (pass/fail)

| Item | Status |
|------|--------|
| \(R^2\) ladder \(0.543 \to 0.746 \to 0.875\) | **PASS** (abstract, table, text, conclusions, best-fit figure) |
| Amp\(=-0.06\), \(\phi\approx0.20\) | **PASS** (ladder text; overlay figure noted as nearby ridge) |
| \(\varepsilon\approx-0.043\), \(\psi\approx176^\circ\) | **PASS** |
| Geometry \(M/M_0=3\), \(d/D=0.5\), \(L\approx0.169\) | **PASS** |
| Ideal integration rescales \(\alpha\) only; expand/collapse flips shape | **PASS** (text + new Fig.) |
| Prefer \(A_{\mathrm{hratio}}\) for sheet comparison | **PASS** |
| Amp\(\times\sin\) / \(2\ell\) labeled diagnostic | **PASS** |
| Nesting dolls speculative | **PASS** |
| CMB adversary present | **PASS** |
| No discovery / cheerleading claim | **PASS** |
| Citations resolve (manual `thebibliography`) | **PASS** (keys match `refs.bib`) |
| Acknowledgments: Rindler + Ackerman; author C. A. | **PASS** |
| Expand/collapse bottom = \(H_{\mathrm{obs}}/H-1\), not B | **PASS** |
| Every `\includegraphics` caption audited vs PNG | **PASS** (with honesty notes above) |

---

## Remaining weaknesses before journal submission

1. **No GR / light-cone / junction construction.** Newtonian skeleton only; Khakshournia-type matching and CMB isotropy with a *beyond-horizon* CoE remain open. A journal referee will demand either a concrete CMB-safe mechanism or a much softer claim.

2. **Phenomenological knobs still carry the high \(R^2\).** Geometry alone is \(R^2\approx0.54\); \(0.875\) needs Amp\(\times\sin\) + \(2\ell\). Consider moving the knobbed panels to an appendix and leading with geometry-only figures.

3. **Illustrative vs best-fit figure parameter drift.** Several figures use \(\kappa=0.231\), \(d/D=0.1\), or \(\phi=0.35\) while the ladder quotes a different best geometry / \(\phi\). Captions now flag this; regenerating a single consistent figure set would be better.

4. **Digitised Migkas errors / covariance not propagated.** RMSE and \(R^2\) treat strip points as independent; strip \(\Delta\ell=90^\circ\) induces strong correlations. Needs a proper likelihood or block-bootstrap before claiming quantitative fit quality.

5. **Bibliography is hand-maintained `thebibliography`.** `refs.bib` exists but is unused by the build; BibTeX/biblatex would reduce drift risk (e.g. Kodama title wording differs slightly between `.bib` and `.tex`).

6. **Related-work depth.** Fennelly / Debney / Pathria entries are thin; journal style may want fuller bibliographic completeness and DOI links.

7. **Tone / venue.** Still a synthesis of an unfinished graduate program. Suitable as a methods / discussion note or arXiv preprint; a high-impact journal will want new data tests or a CMB-compatible metric, not only the Newtonian \(R^2\) ladder.

8. **Unused asset.** `figures/expand_vs_collapse_vs_migkas.png` remains in the tree (invsq + Migkas B). Safe to delete or quarantine so it is not re-included by mistake.

---

## Top 5 editorial changes (summary)

1. Replaced expand/collapse figure; bottom panel correctly \(H_{\mathrm{obs}}/H_{\mathrm{bg}}-1\).
2. Honest A/B caption on Migkas rotation overlay.
3. Clarified AF / rotation captions so illustrative slices ≠ ladder peak.
4. Tightened abstract + conclusions; removed informal overclaim language.
5. Rebuilt `main.pdf` (14 pp); wrote this review.

