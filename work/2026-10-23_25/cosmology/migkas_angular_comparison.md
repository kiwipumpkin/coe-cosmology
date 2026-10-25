# Migkas angular match — problem list

**Date:** 2026-09-18 (CT)  
**Success criterion (Chris):** first-principles \(H_{\mathrm{obs}}(\hat n)\) (mapped to \(A(\ell)\)) matches Migkas **angular / sky** data — **not** a best-fit \(A\cos\varphi\) dipole.

**Sources:** Migkas et al. 2020, A&A 636, A15 (arXiv:2004.03305). Digitised from rendered Fig. 6 (A top, B bottom vs Galactic longitude, \(\Delta\ell=90^\circ\) bins) and Fig. 7 top-left (Our sample, B fixed).

---

## What Migkas publishes (angular pattern)

| Item | Result |
|------|--------|
| Method | Sliding sky cones / longitude bins; fit \(L_X\!-\!T\) normalisation \(A\) (or \(H_0\) with \(A\) fixed) |
| 1D plot (Fig. 6) | \(A/A_{\mathrm{all}}\) and \(B/B_{\mathrm{all}}\) vs **Galactic longitude** (\(\Delta\ell=90^\circ\), full \(b\)) |
| Low-\(A\) / high-\(H_0\) | Their sample: toward \((l,b)\sim(280^\circ,-20^\circ)\); joint 3-sample: \((303^\circ,-27^\circ)\) |
| Amplitude | Max \(D_L\) contrast \(\sim 16\pm 3\%\) (joint); \(A\) varies by tens of percent in 1D |
| Shape (their words / maps) | Strong directional anisotropy; often summarised with a **dipole** scan, but 1D \(A(\ell)\) has one deep trough near \(\ell\sim 250^\circ\)–\(300^\circ\) and a broad high near \(\ell\sim 0^\circ\)–\(60^\circ\) (more **monopole+dipole-like** than pure \(P_2\)) |
| Slope \(B\) | Much closer to isotropic than \(A\) (Fig. 6 bottom) |

**Conversion used here:** \(A\times H_0^2\) conserved \(\Rightarrow\) \(H/H_{\mathrm{ref}}=1/\sqrt{A/A_{\mathrm{all}}}\).

Digitised CSVs: `migkas/digitised/migkas2020_fig6_{A,B}_vs_longitude.csv`, `migkas2020_fig7_OurSample_A_vs_longitude.csv`.

---

## What our model predicts (binding)

| Item | Prediction |
|------|------------|
| Ontology | CoE beyond horizon; observer at \(D\); decreasing \(\rho(r)\); Gauss \(g=-GM(<r)/r^2\) |
| Homology (\(M=0\)) | \(H_{\mathrm{obs}}=H\) isotropic \(\Rightarrow A/A_{\mathrm{all}}=1\) |
| Gravity, \(d\ll D\) | Leading **quadrupole** \(P_2(\mu)\), \(\mu=\hat n\cdot\hat D\); \(\kappa_{\mathrm{nat}}=0.15\) |
| Full-sky multipoles | \(C_2\sim 2\alpha\kappa=0.30\) (bolt-on); \(C_1\approx 0\); small \(C_3\propto d/D\) |
| Mapped to \(A(\ell)\) | Prefer axis at Migkas low-\(A\) longitude \(\ell\approx 280^\circ\); \(A/A_{\mathrm{all}}=1/(1+f)^2\) renormalised to mean 1 |

**Figure:** `figures/Hobs_vs_angle.png` — model (raw + \(\Delta\ell=90^\circ\)-smoothed) + digitised Migkas Fig. 6; **no** \(A\cos\varphi\) overlay.

---

## Concrete problems / mismatches to fix (ordered attack list)

1. **Morphology: \(P_2\) vs Migkas nearly-dipole \(A(\ell)\)**  
   Model (even after \(90^\circ\) smoothing) has **two** low-\(A\) troughs (toward \(\pm\hat D\)). Migkas Fig. 6 has **one** deep low-\(A\) trough near \(\sim 280^\circ\) and one high region near \(\sim 30^\circ\).  
   **Fix direction:** odd multipoles / bulk flow / depth \(d/D\) large enough for \(C_1\) or \(C_3\); selection weighting; or cone response of \(P_2\) through Migkas exact pipeline (not longitude strip only).

2. **Observational filter: longitude strip ≠ their 2D cone maps**  
   Fig. 6 averages over all \(b\) in \(\Delta\ell=90^\circ\). True Migkas signal is in **2D cone scans** (Figs. 8–9, 21–23).  
   **Fix:** Forward-model \(H_{\mathrm{obs}}(\hat n)\) through **identical \(\theta=75^\circ\) cone scan** and compare to their \(H_0\)/ \(A\) maps (digitise map or reimplement scan on public cluster table).

3. **Amplitude after fair smoothing still not tuned to data**  
   At \(\kappa_{\mathrm{nat}}=0.15\), \(\alpha=1\), unsmoothed \(A\) swings are large; \(\Delta\ell=90^\circ\) smoothing reduces them but may still overshoot digitised Fig. 6. Drag \(\alpha\sim 0.66\) and/or enclosed-mass \(\kappa<\kappa_{\mathrm{nat}}\) need a **pipeline-matched** amplitude check (not dipole MCMC).

4. **Axis / latitude: \(b\neq 0\)**  
   Preferred direction is at \(b\sim -20^\circ\) to \(-27^\circ\), not on the plane. Longitude-only model cut is schematic.  
   **Fix:** Evaluate model on full sky with \(\hat a=(303^\circ,-27^\circ)\) (joint) or \((280^\circ,-20^\circ)\) (their sample) and project with their weights.

5. **\(B(\ell)\) is a null test we currently ignore**  
   Migkas \(B\) is nearly flat; our model does not predict \(L_X\!-\!T\) slope anisotropy. Keep \(B\) as a **must-not-break** diagnostic when adding selection / ICM effects.

6. **Digitisation honesty**  
   Points are from pixel centroids on the arXiv page render (±~0.02–0.03 systematic). Prefer author tables / VizieR cluster list + re-fit for publication-grade comparison.

---

## Ordered attack list (Chris goal)

1. Implement **Migkas-like \(\theta=75^\circ\) cone forward map** of model \(H_{\mathrm{obs}}(\hat n)\) (and \(A\propto H^{-2}\)) with axis at \((280^\circ,-20^\circ)\) and \((303^\circ,-27^\circ)\).  
2. Quantify **shape residual**: dipole vs quadrupole power in model-through-pipeline vs digitised / re-fit Migkas maps (problem #1).  
3. Add **odd terms** only if motivated (finite \(d/D\), bulk flow, or asymmetric selection) — not by fitting \(A\cos\varphi\) by hand.  
4. Retune **\(\alpha,\kappa\)** only after pipeline match; compare to \(\kappa_{\mathrm{nat}}=0.15\).  
5. Keep **\(B\) flat** and galactic-absorption tests as vetoes.

**Next fix to implement:** (1) cone-scan forward model matching Migkas §4.3 / Fig. 21–23.

---

## Drive deliverables (this turn)

- Revised `Hobs_vs_angle` PNG (Migkas Fig.6 overlay, no dipole)
- Combined PDF (revised)
- Digitised CSVs (Fig.6 A,B; Fig.7 Our sample A)
- This note + Google Doc twin
