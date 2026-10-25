# Migkas angular overlay — deliverable report

**Time:** 2026-09-18 ~10:20 CT  
**Goal:** First-principles model vs Migkas **angular** data (Fig.6 A,B vs Galactic longitude). **No** \(A\cos\varphi\) dipole overlay.

## Done

1. **Digitised** Migkas+2020 Fig.6 A/Aall and B/Ball vs ℓ from arXiv PDF page render (`/workspace/cosmology/migkas/migkas2020.pdf`).
2. **Digitised** Fig.7 top-left (Our sample, B fixed) A/Aall vs ℓ.
3. **Updated** `figures/make_figures.py` `fig_Hobs_vs_angle`: homology + ST+Gauss model (raw + Δℓ=90° smoothed) mapped to A/Aall with preferred axis at ℓ=280°; Migkas points overlaid; **no dipole template**.
4. Regenerated PNG + combined PDF.
5. Wrote problem list markdown + Google Doc.
6. **Uploaded** revised PNG, PDF, schematic, and DIGITISED CSVs to Drive (folder below).

## Drive links (folder `1wZYStJhpdowbcjTkkS9TOImAtt9_GvzD`)

| Item | Link |
|------|------|
| Problem list (Google Doc) | https://docs.google.com/document/d/1XxbxwWEAPT0DYBfyd9FOZ-m7hu2MeR1jdv3mXDxJvC0/edit |
| **Hobs vs angle PNG** (revised, Migkas overlay) | https://drive.google.com/file/d/14rGW_aTnKSSEl2ZDP_5mF-v4Yj53-bYe/view |
| **Combined figures PDF** (revised) | https://drive.google.com/file/d/1RnidJS3PrKv99TCExFHxMqy_lTyWB7w6/view |
| Model vs Migkas schematic PNG | https://drive.google.com/file/d/1l-R8fdJ4eJC7wbReuxoxmMApkWjjj8fG/view |
| Fig.6 A CSV (DIGITISED, re-upload) | https://drive.google.com/file/d/14Jh75HYrKiTseFRLwUmy5Q7-kfEhf05i/view |
| Fig.6 B CSV (DIGITISED, re-upload) | https://drive.google.com/file/d/1smgsFks0NckcTekJinB0Wh0u9zOmZTJH/view |
| Fig.7 OurSample A CSV (DIGITISED, re-upload) | https://drive.google.com/file/d/1e1BZlSD5uunv3wLnhwsJJlbzRQpTpqUI/view |

Earlier CSV uploads (same content, older IDs):  
A `1xmXtIjQF6jeWoQaedjkehOS875IKNbN5` · B `1w854teiiCHdVaaVRBtoe0UPO4WTBERtq` · Fig7 `1sc-YnKgV-QsHUui9Wh2xmRQxfN1n_3Hf`

**Local staging:** `/workspace/cosmology/figures/revised_no_dipole/` and `/workspace/cosmology/migkas/digitised/`

## Top 5 problems to fix

1. **Morphology:** Model \(P_2\) → two low-A troughs; Migkas Fig.6 → one deep trough ~250–300° (dipole-like).
2. **Filter:** Longitude strip ≠ Migkas 2D cone maps — need θ≈75° cone forward model.
3. **Amplitude:** After Δℓ=90° smoothing, retune α,κ only vs pipeline, not dipole MCMC.
4. **Axis latitude:** Prefer â=(280°,−20°) or (303°,−27°), not b=0 cut.
5. **B null:** Keep Migkas near-flat B as must-not-break diagnostic.

## Next fix to implement

**Migkas-like θ=75° cone-scan forward map** of model \(H_{\mathrm{obs}}(\hat n)\) / \(A\propto H^{-2}\) with axis at (280°,−20°) and (303°,−27°).
