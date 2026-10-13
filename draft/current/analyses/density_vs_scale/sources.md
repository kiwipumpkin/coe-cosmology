# Density contrast vs region size: points and sources

δ = ρ/ρ_mean − 1 (matter, unless noted). Radii are converted to Mpc with h = 0.70. "Derived" means I computed δ from the paper's mass and volume using ρ_mean = Ω_m ρ_crit with Ω_m = 0.3 (≈ 8.3×10¹⁰ h² M☉/Mpc³). Those values are approximate (±30% or worse).
Chart: `density_vs_scale.png`. Script: `make_chart.py`. ΛCDM σ(R) table: `lcdm_sigma_tophat.txt`.

## Plotted, over-dense
| point | radius (Mpc) | δ | how | source |
|-|-|-|-|-|
| Coma cluster, R200 | 2.84 (1.99 h⁻¹) | ≈ +645 | Weak lensing M200 = 1.88×10¹⁵ h⁻¹ M☉. By definition, R200 encloses 200× the critical density. I assumed critical, not mean, which gives 200/Ω_m − 1 ≈ 645 (approx.). Cluster cores (R2500) run ~10× higher still. | Kubo et al. 2007, ApJ 671, 1466, doi:10.1086/523101 |
| Shapley core | 11.4 (8 h⁻¹) | ≈ +20 (derived). Range of earlier estimates +9 to +66 | 3.8×10¹⁵ h⁻¹ M☉ inside an 8 h⁻¹ Mpc sphere. Reisenegger 2000 gave ρ/ρ_c ≈ 3–20. | Ragone et al. 2006, A&A 445, 819, doi:10.1051/0004-6361:20053623 |
| Saraswati supercluster | 36.6 (radius of a sphere with the convex-hull volume 2.05×10⁵ Mpc³) | ≳ +1.62 (lower limit, quoted by the authors) | M ≈ 2×10¹⁶ M☉ in 43 clusters, at z ≈ 0.28 | Bagchi et al. 2017, ApJ 844, 25, arXiv:1707.03082 |
| Sloan Great Wall | 79 (equivalent sphere for V = 7.2×10⁵ h⁻³ Mpc³) | ≈ +1.0 (derived) | M ≈ 1.2×10¹⁷ h⁻¹ M☉ (±20%). The authors say it is hard to reconcile with Gaussian ΛCDM if σ8 < 0.9. | Sheth & Diaferio 2011, MNRAS 417, 2938, arXiv:1105.3378 |

## Plotted, under-dense
| point | radius (Mpc) | δ | how | source |
|-|-|-|-|-|
| Boötes void | ≈ 44 (radius ≈ 3,100 km/s; the paper quotes "62 Mpc" for its older H0) | ≈ −0.67 to −0.83 (galaxy density 1/3 to 1/6 of mean) | Galaxy counts, not mass. The conversion of the radius is approximate. | Kirshner et al. 1987, ApJ 314, 493, doi:10.1086/165080; density range as quoted in Szomoru et al. 1996, AJ 111, 2150 |
| SDSS voids, stacked lensing | 21–79 (15–55 h⁻¹) | ≈ −0.4 inside R_v | Weak lensing of about 20,000 voids. Mass, not galaxies. Typical voids, not the extreme ones. | Clampitt & Jain 2015, MNRAS 454, 3357, arXiv:1404.1834 |
| Eridanus supervoid (Cold Spot) | 314 (220 ± 50 h⁻¹) | −0.14 ± 0.04 (matter) | Galaxy δ ≈ −0.2 divided by bias. The authors call it a ≥3.3σ ΛCDM fluctuation. Kovács & García-Bellido 2016 (arXiv:1511.09008) give δ ≈ −0.2 to −0.25 at R ≈ 195 h⁻¹. | Szapudi et al. 2015, MNRAS 450, 288, arXiv:1405.1566 |
| Local Hole / KBC (open circle, **disputed**) | ≈ 215–300 | ≈ −0.15 to −0.33 (galaxies/light) | Keenan et al.: local K-band light density about 1.5× lower than at z > 0.1, out to ~300 Mpc, which implies ≈ −0.3. Whitbourn & Shanks: −15 ± 3% to 150 h⁻¹ Mpc. **Rebuttal:** Kenworthy et al. 2019 find SN Ia distances exclude \|δ\| > 20% local voids at 4–5σ. | Keenan, Barger & Cowie 2013, ApJ 775, 62, arXiv:1304.2884; Whitbourn & Shanks 2014, MNRAS 437, 2146, arXiv:1307.4405; Kenworthy, Scolnic & Riess 2019, ApJ 875, 145, arXiv:1901.08681 |

## Plotted as context (not a density)
- **Quasar count dipole** (open diamond): amplitude D = 0.0155 from 1.36 M CatWISE quasars, against ~0.007 expected from our motion. That is roughly 3% end to end across the sky, 4.9σ. It is a tilt in counts, not a measured δ. The x position (~7 Gpc) is only schematic, because the sample spans z ≈ 0.1–3. Source: Secrest et al. 2021, ApJL 908, L51, arXiv:2009.14826.
- **Homogeneity scale** (dotted line at ~95 Mpc): the fractal dimension reaches 3 within 1% above 70 ± 5 h⁻¹ Mpc (WiggleZ) and 63.3 ± 0.7 h⁻¹ Mpc (BOSS; the ΛCDM prediction is 62.0). Sources: Scrimgeour et al. 2012, MNRAS 425, 116, arXiv:1205.6812; Ntelis et al. 2017, JCAP 06, 019, arXiv:1702.02159.
- **Size-only giant structures** (grey arrow, ~200–1500 Mpc). They are claimed as clusterings of tracers, and none has a published mass δ:
  - Giant Arc: ~1 Gpc long, z ≈ 0.8, 3.0–4.8σ. Lopez et al. 2022, MNRAS 516, 1557, arXiv:2201.06875.
  - Big Ring: ~400 Mpc diameter, up to 5.2σ. Lopez et al. 2024, JCAP 07, 055, arXiv:2402.07591.
  - Huge-LQG: ~500 Mpc characteristic, longest ~1.24 Gpc, z ≈ 1.27. Clowes et al. 2013, MNRAS 429, 2910, arXiv:1211.6256.
  - Hercules–Corona Borealis Great Wall: GRB clustering, ~2–3 Gpc, z ≈ 1.6–2.1. Horváth et al. 2015, A&A 584, A48, arXiv:1510.01933.

## Not plotted (checked, left out for clarity)
- **Shapley whole**: 1.8×10¹⁶ h⁻¹ M☉ in 1.2×10⁵ h⁻³ Mpc³, which gives δ ≈ +0.8 at R ≈ 44 Mpc (derived; Sheth & Diaferio 2011). It sits right on Saraswati.
- **BOSS Great Wall** (z ≈ 0.47): superclusters 64–186 h⁻¹ Mpc across. The luminosity density is ~9–10× mean, but only in 8 h⁻¹ Mpc smoothing cells, not averaged over the whole wall. Total mass 4.2×10¹⁶ h⁻¹ M☉. Lietzen et al. 2016, A&A 588, L4, arXiv:1602.08498; Einasto et al. 2017, arXiv:1703.08444.
- **Laniakea**: ~160 Mpc across, ~10¹⁷ M☉. It is defined by flow boundaries, not by density, so its δ is not well defined. A rough derived value is only ~+0.1. Tully et al. 2014, Nature 513, 71, arXiv:1409.0880.
- **Void cores in simulations**: universal profile central δ_c ≈ −0.95 at z = 0. These are simulations, not observations. Hamaus, Sutter & Wandelt 2014, PRL 112, 251302, arXiv:1403.5499.

## Model markers (from our own fits, not observations)
- **Required smooth field** (`../inside_visible/smooth_field/summary.md`), Gaussian smoothing σ_s:
  - 1 Gpc: peak 0.068, rms 0.0105
  - 2 Gpc: peak 0.055, rms 0.012
  - 4 Gpc: peak 0.059, rms 0.014
  - These are plotted at x = σ_s. A Gaussian of width σ_s averages over more volume than a top-hat of radius σ_s, so ΛCDM's expectation for these points is really ~5× *lower* than the band drawn. The chart is generous to the field.
- **Point-mass fit cores**: δ = +100, radius ≈ 0.83 Gpc for a 2×10⁵² kg core (from M = 100·ρ̄·V with ρ̄ = 2.9×10⁻²⁷ kg/m³).
- **Point-mass fit voids**: δ = −0.7, R_v 3.5–4 Gpc (`../inside_visible/summary.md`).

## ΛCDM band
Linear σ(R) for a top-hat sphere at z = 0: Eisenstein & Hu 1998 no-wiggle transfer function, σ8 = 0.81, Ω_m = 0.31, Ω_b = 0.049, h = 0.70, n_s = 0.965. Code: `../inside_visible/smooth_field/lcdm.py`.

| R (Mpc) | 10 | 100 | 300 | 1000 | 2000 | 4000 |
|-|-|-|-|-|-|-|
| σ | 0.89 | 0.10 | 0.021 | 0.0028 | 0.00079 | 0.00022 |

Below ~10 Mpc, σ > 1 and linear theory badly underestimates real (nonlinear) extremes, so that region is greyed out. The ±3σ band for under-densities is not Gaussian-accurate near δ = −1.

---
# v2 additions (Oct 9, 2026): more points, used in `density_vs_scale_fits_radius_v2.png`, `fits_v2.py` and `points_v2.json`
All radii are in Mpc with h = 0.7. "Derived" means δ = M/(ρ̄V) − 1 with Ω_m = 0.3. Tracer is given in brackets.

## Over-dense, new
| point | R (Mpc) | δ | tracer / how | source |
|-|-|-|-|-|
| Perseus R200 (definitional) | 1.79 | ≈ +645 (200ρ_c assumed) | X-ray (Suzaku) r200 = 1.79 ± 0.04 Mpc | Simionescu et al. 2011, Science 331, 1576, arXiv:1102.2429 |
| Virgo zero-velocity sphere | 7.2 ± 0.7 | ≈ +11.5 (derived; +8 to +15 from the mass error) | mass from infall: M = (8.0 ± 2.3)×10¹⁴ M☉ | Karachentsev et al. 2014, ApJ 782, 4, arXiv:1312.6769 |
| Corona Borealis SC | 17.9 (turnaround radius 12.5 h⁻¹) | ≈ +37, range +8 to +175 (derived from M = 0.6–12×10¹⁶ h⁻¹ M☉) | dynamics + caustics, all from one paper | Pearson, Batiste & Batuski 2014, MNRAS 441, 1601, arXiv:1404.1308 |
| BOSS Great Wall A / B / C / D | 36.1 / 36.7 / 23.7 / 21.2 (equivalent spheres for V = 67,500 / 70,848 / 19,008 / 13,635 h⁻³ Mpc³) | +2.7 / +0.87 / +3.4 / +1.6 (derived) | masses 2.1 / 1.1 / 0.7 / 0.3 ×10¹⁶ h⁻¹ M☉, from luminosity × M/L | Einasto et al. 2017, A&A 603, A5, arXiv:1703.08444 |
| Shapley whole (moved onto the chart) | 43.7 | +0.80 (derived) | mass | Sheth & Diaferio 2011, arXiv:1105.3378 |
| Horologium–Reticulum SC | ≈ 35 (rough: 12°×12° and 17,000–22,500 km/s) | ≈ +1.4 (inter-cluster regions only) | galaxy counts | Fleenor et al. 2005, AJ 130, 957, arXiv:astro-ph/0505361 |
| Vela SC | ≈ 50 (rough: extent ~115×90 Mpc) | +0.50 to +0.77 (plotted +0.64) | galaxy counts (2MASS photo-z shell); preliminary | Kraan-Korteweg et al. 2017, MNRAS 466, L29, arXiv:1611.04615 |
| **2M++ densest sphere** (my calculation from the public map) | 14.3 / 28.6 / 42.9 / 71.4 / 100 / 143 (10–100 h⁻¹) | +6.3 / +1.9 / +1.0 / +0.35 / +0.21 / +0.082 (mass) | see the 2M++ note below | Carrick et al. 2015, MNRAS 450, 317, arXiv:1504.04627; data at cosmicflows.iap.fr |

## Under-dense, new
| point | R (Mpc) | δ | tracer / how | source |
|-|-|-|-|-|
| Local Void | ≈ 23 (Tully: void centre ≥ 23 Mpc away) | ≈ −0.8 (galaxy number density inside is ~1/5 of outside) | galaxies | Nasonova & Karachentsev 2011, Astrophysics 54, 1, arXiv:1011.5985; Tully et al. 2008, ApJ 676, 184, arXiv:0705.4139 |
| BOSS DR12 largest void | 91.5 (R_eff 64.05 h⁻¹) | −0.89 (**central minimum**, not the average inside R) | galaxies (ZOBOV); catalog-wide mean central density is ~0.3× mean | Mao et al. 2017, ApJ 835, 161, arXiv:1602.02771 |
| Eridanus, alternative fit | 279 (transverse 195 h⁻¹; line of sight ~500 h⁻¹) | −0.25 (**central** δ₀) | matter (galaxies ÷ bias) | Kovács & García-Bellido 2016, MNRAS 462, 1882, arXiv:1511.09008 |
| Local Hole: Böhringer 2020 (disputed family) | 100 / 140 | −0.30 ± 0.15 / −0.20 ± 0.10 | X-ray clusters, bias-corrected to matter; they quote a 1.3–3.7σ ΛCDM tension | Böhringer, Chon & Collins 2020, A&A 633, A19, arXiv:1907.12402 |
| Local Hole: Wong 2022 (disputed) | ≈ 320 (z < 0.075) | −0.21 (2M++/2MRS n(z)); −0.20 ± 0.02 from n(K) | galaxies, K-band | Wong et al. 2022, MNRAS 511, 5742, arXiv:2107.08505 |
| Local Hole: Keenan 2013 and Whitbourn & Shanks 2014 | 300 / 214 | −0.33 / −0.15 | now plotted as two separate points (v1 had one combined point) | see above |
| **2M++ emptiest sphere** (my calculation) | 14.3 / 28.6 / 42.9 / 71.4 / 100 / 143 | −0.75 / −0.59 / −0.45 / −0.23 / −0.13 / −0.041 (mass) | see the 2M++ note below | Carrick et al. 2015 |

**2M++ note.** I downloaded the public 2M++ density cube (257³, Galactic coordinates ±200 h⁻¹ Mpc, luminosity-weighted δ_g*, already smoothed with a 4 h⁻¹ Mpc Gaussian). Then:
- I averaged it in top-hat spheres of radius R, keeping only centres with r + R ≤ 150 h⁻¹ Mpc.
- I converted to mass with δ_m = δ_g*/b*, where b* = Ω_m^0.55/β* = 1.20 (β* = 0.43).
- The densest sphere at R ≤ 30 h⁻¹ is centred on Shapley, as a sanity check.
- Caveat 1: at large R only centres near us qualify. At 100 h⁻¹ that means centres within 50 h⁻¹ of us, so the largest-R extremes come from a small volume and are biased low.
- Caveat 2: the nested spheres are strongly correlated, so the 2M++ fit errors are underestimated.
- Script: `twompp_extremes.py`. Output: `twompp_extremes.json`. Data: `data/twompp_density.npy`.
- 2M++ rms of mass δ: 0.48 / 0.31 / 0.22 / 0.12 / 0.062 / 0.022 at 10–100 h⁻¹. Compare ΛCDM linear σ = 0.10 at 100 Mpc.

## Checked, not plotted (no usable δ with a size)
- **Dipole Repeller** (Hoffman et al. 2017, Nature Astronomy 1, 0036, arXiv:1702.02483): inferred from velocities only, no δ quoted.
- **Hercules SC**: no clean δ found.
- **Einasto/Liivamägi supercluster catalogs**: they use luminosity-density thresholds (D8 ≥ 5 in 8 h⁻¹ Mpc cells), not whole-object averages. The BOSS Great Wall D8 means of 9–10 fall in this class.
- **Hamaus 2014/2016 void profiles**: the 2014 values come from simulations; the 2016 BOSS profiles are fit parameters by size bin, with no single extreme.
- **Laniakea**: defined by flows, so no δ.
- **Count-in-cells extremes** (Gaztañaga etc.): no published maximum-cell values at Mpc–Gpc scales found quickly. The 2M++ calculation above stands in for them.
- **Homogeneity scale** (Scrimgeour 2012; Ntelis 2017; Gonçalves 2021 for quasars): these give a scale where counts are within 1% of uniform, not an extreme δ. It stays as context only.
