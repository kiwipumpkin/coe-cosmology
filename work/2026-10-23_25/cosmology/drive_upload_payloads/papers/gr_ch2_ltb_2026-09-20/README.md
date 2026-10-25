# GR Chapter 2 — Off-center LTB sketch (2026-09-20)

Standalone LaTeX note promoting the Newtonian ST+CoE Migkas sky template to dust LTB (optional Λ) language. First shot, not a finished GR cosmology.

## File map

| Path | Role |
|------|------|
| `main.tex` | Chapter note |
| `refs.bib` | Bibliography |
| `main.pdf` | Compiled PDF |
| `figures/A_vs_ell_toy_ltb.png` | Toy \(A(\ell)\) vs digitised Migkas Fig. 6 |
| `figures/cmb_dipole_vs_offset.png` | CMB dipole vs observer offset (kill-shot sketch) |
| `figures/toy_density_H_profile.png` | Toy \(\rho(r)\) and FLRW-shell \(H(r)\) |
| `scripts/toy_offcenter_ltb.py` | Numerical toy (matplotlib) |
| `toy_results.json` | Numbers from last toy run |
| `NOTES_FOR_CHRIS.md` | Plain-speech summary |
| `README.md` | This file |

## Rebuild

```bash
cd /workspace/cosmology/papers/gr_ch2_ltb_2026-09-20

# regenerate figures
/workspace/cosmology/.venv/bin/python scripts/toy_offcenter_ltb.py

# compile PDF
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Migkas input CSV (do not invent numbers):

`/workspace/cosmology/migkas/fig6_digitized.csv`

## Honesty labels

- Toy uses FLRW-shell \(H(r)\) from enclosed \(M(r)\), not full LTB null geodesics.
- Sheet AF: \(A\propto H_\mathrm{obs}/H\) (Chapter 1 convention).
- Amp×sin / 2ℓ are diagnostic residuals, not GR parameters.
- CMB curve is linear \(|\nabla\ln H|\,d\) plus Alnes-scaled reference; ceiling class ~15 Mpc (void) or tighter for this steep profile.
