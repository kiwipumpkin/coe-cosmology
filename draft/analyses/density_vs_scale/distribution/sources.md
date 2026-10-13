# Sources: density distribution vs region size (Task 4)

## Data used
| Data | What it gives | Reference |
|---|---|---|
| 2M++ density field (public 257³ cube, ±200 h⁻¹ Mpc) | δ for every sphere of radius 10–150 Mpc (counts-in-cells) | Carrick et al. 2015, MNRAS 450, 317, arXiv:1504.04627 |
| BOSS DR12 void catalog (VizieR J/ApJ/835/161), 1,228 voids | void radius Reff, central δ (delmin) | Mao et al. 2017, ApJ 835, 161, arXiv:1602.02771, doi:10.26093/cds/vizier.18350161 |
| BOSS LOWZ/CMASS footprint areas | survey volume for void counts | Reid et al. 2016, MNRAS 455, 1553, arXiv:1509.06529 |
| redMaPPer SDSS DR8 clusters (VizieR J/ApJ/785/104), λ ≥ 20 | cluster richness and redshift | Rykoff et al. 2014, ApJ 785, 104, arXiv:1303.3562, doi:10.1088/0004-637X/785/2/104 |
| Richness to mass: M200m = 10^14.344 (λ/40)^1.33 h⁻¹ M☉ | cluster radius R200m | Simet et al. 2017, MNRAS 466, 3103, arXiv:1603.06953 |
| SDSS DR7 superclusters (VizieR J/A+A/539/A80, main50), 982 objects | supercluster volume → equivalent radius | Liivamägi, Tempel & Saar 2012, A&A 539, A80, arXiv:1012.1989, doi:10.1051/0004-6361/201016288 |
| Halo mass function (model) | theory cluster counts | Tinker et al. 2008, ApJ 688, 709, arXiv:0803.2706 |
| Void size function, Vdn (model) | theory void counts | Jennings, Li & Hu 2013, MNRAS 434, 2167, arXiv:1304.6087 |
| Linear ΛCDM σ(R) | standard-cosmology width | Eisenstein & Hu 1998 no-wiggle P(k), σ8 = 0.81, Ωm = 0.31, h = 0.70 (our lcdm.py) |
| Extreme points and fits | the 37 points and fits from v2 | ../sources.md, ../fits_v2.json |

## Derived or approximate (marked as such)
- **2M++ mass δ = galaxy δ / 1.20 (DERIVED).** The bias b* = Ωm^0.55/β* = 0.3^0.55/0.43.
  - The cube already carries 4 h⁻¹ Mpc Gaussian smoothing, so the widths at R ≤ 14 Mpc are too narrow.
  - Only centres with r + R ≤ 150 h⁻¹ Mpc are used.
  - The number of independent spheres is about (150 h⁻¹/R)³: ~9,800 at 10 Mpc, ~19 at 80 Mpc and ~3 at 150 Mpc. The 150 Mpc histogram is really about one region (ours).
  - The local volume is over-dense, so the mean δ is +0.04 at 80–150 Mpc.
  - 2M++ is a reconstruction that shrinks toward the mean at large distance. Its widths at R ≥ 110 Mpc are therefore biased low.
- **Volumes (APPROXIMATE).** These use the flat Ωm = 0.31 comoving distance and the footprint areas, with no masks or completeness corrections.
  - redMaPPer area: 10,400 deg² (approximate).
  - Mao voids are the *cleaned* sample, so the counts are lower limits.
- **Radii (DERIVED).** Void R = Reff/h. Cluster R = R200m from the Simet relation; individual masses carry about 0.25 dex scatter. Supercluster R = (3 Vol/4π)^(1/3)/h.
- **ΛCDM histogram curve (MODEL).** Lognormal with the linear σ. The σ is convolved with the same 4 h⁻¹ Gaussian so it matches 2M++. Nonlinear effects at R ≤ 20 Mpc are only approximated.
- **Width trend (DERIVED FIT).** The 2M++ std is fitted at R = 20–80 Mpc and extrapolated to 1–4 Gpc. This is far outside the data.
