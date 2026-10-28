# CoE data fetch manifest — 2026-10-27 (CDT)
All files in /workspace/coepdf/data_fetch/. Sizes in bytes; sha256 prefix shown.

## 1. Pantheon+ — GOT
Source: https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/
- pantheon/Pantheon+SH0ES.dat — 579283 B, 1701 SNe rows. Columns: CID IDSURVEY zHD zHDERR zCMB zCMBERR zHEL zHELERR m_b_corr m_b_corr_err_DIAG MU_SH0ES MU_SH0ES_ERR_DIAG CEPH_DIST IS_CALIBRATOR USED_IN_SH0ES_HF c cERR x1 x1ERR mB mBERR x0 x0ERR COV_* RA DEC HOST_RA HOST_DEC HOST_ANGSEP VPEC VPECERR MWEBV HOST_LOGMASS ... (full header in file)
- pantheon/Pantheon+SH0ES_STAT+SYS.cov — 33284960 B (33 MB, not 60 MB). First line = 1701 (N), then N² entries.

## 2. Planck PR3
Source dir: https://irsa.ipac.caltech.edu/data/Planck/release_3/
- COM_CMB_IQU-smica_1024_R3.00_full.fits — **FAILED**: IRSA .../maps/component-maps/cmb/COM_CMB_IQU-smica_1024_R3.00_full.fits → HTTP 404 (not in IRSA listing); PLA http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=COM_CMB_IQU-smica_1024_R3.00_full.fits → HTTP 500.
- planck/COM_CMB_IQU-smica_2048_R3.00_full.fits — GOT (fallback), 2013312960 B (1.9 GB). Header OK: NSIDE 2048, NESTED, GALACTIC; cols I_STOKES Q_STOKES U_STOKES TMASK PMASK I_STOKES_INP Q_STOKES_INP U_STOKES_INP TMASKINP PMASKINP (K_CMB).
- planck/smica_I_nside1024_udgrade_from2048.fits — DERIVED here (not an official product): I_STOKES ud_graded to Nside 1024, RING, 50339520 B. Script: planck/make_smica1024.py.
- planck/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits — GOT, 201335040 B, from .../ancillary-data/masks/. NSIDE 2048, NESTED, GALACTIC, col TMASK.
- planck/COM_PowerSpect_CMB-TT-full_R3.01.txt — copied from letter folder, 170547 B.

## 3. Clusters
### eROSITA eRASS1 (Bulbul et al. 2024, A&A 685, A106) — GOT
Source: https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/ (.fits direct → HTTP 403; the .fits.tgz versions downloaded and extracted). The path in the request (.../Catalogues/) → HTTP 404.
- erosita/erass1cl_primary_v3.2.fits — 50636160 B, 12,247 clusters. Key cols: RA DEC BEST_Z BEST_ZERR L500 L500_L L500_H (0.2–2.3 keV) Lbol500 KT KT_L KT_H YX500 M500 R500 … (87 cols).
- erosita/erass1cl_cosmology_v1.1.fits — 1186560 B, 5,259 clusters (cosmology subsample, Ghirardini+2024). Cols: RA DEC Z_LAMBDA Z_LAMBDA_E KT KT_L KT_H L500 L500_L L500_H CR500 … MGAS500 YX500.
- CSV exports (*_lbz.csv): NAME RA DEC GLON GLAT z L500 L500_L L500_H KT KT_L KT_H YX500 YX500_L YX500_H (GLON/GLAT via astropy ICRS→Galactic; script erosita/make_csv.py). _L/_H are lower/upper bounds as given in the catalogue, not ±errors.
- CAVEAT: KT = -1 is a no-measurement flag. Primary: only 3,389 of 12,247 have KT>0. Cosmology: 5,023 of 5,259 have KT>0. No Y_SZ column (YX500 is X-ray Y_X, not SZ). Units per data-model page erass1cl_primary_v3.2.html (L500 in 1e42 erg/s per that page — check before use).
### Migkas et al. 2021 SZ sample (A&A 649, A151) — FAILED
- http://cdsarc.cds.unistra.fr/ftp/J/A+A/649/A151/ReadMe → HTTP 500; http://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/A%2BA/649/A151 → HTTP 500; https to cdsarc/vizier (incl. IUCAA mirror) → TLS error (unexpected EOF); CfA mirror → timeout. CDS is blocked from this box too. (VizieR catalogue ID J/A+A/649/A151 not verified because CDS was unreachable.)

## 4. Number-count dipoles (from the papers, read 2026-10-27)
- Secrest et al. 2021, ApJL 908, L51 (arXiv:2009.14826), CatWISE2020 quasars, 1,355,352 objects: **D = 0.01554 toward (l,b) = (238.2°, +28.8°)**, 27.8° from CMB dipole; kinematic expectation ≈ 0.007; p = 5e-7 (4.9σ).
- Secrest et al. 2022, ApJL 937, L31 (arXiv:2206.05624):
  - NVSS (S>10 mJy, 508,144 sources): **D = (1.23 ± 0.25)e-2 toward (l,b) = (196° ± 13°, +46° ± 10°)**, 45° from CMB dipole; kinematic expectation 0.41e-2; 2.6σ.
  - WISE quasars (W1<16.5): **D = (1.48 ± 0.16)e-2 toward (238° ± 7°, +31° ± 5°)**, 26° off; expectation 0.73e-2; 4.4σ.
  - Joint common dipole: (1.40 ± 0.13)e-2 toward (233° ± 6°, +34° ± 5°); after removing kinematic part: (0.86 ± 0.14)e-2 toward (217° ± 10°, +20° ± 7°). Joint 5.1σ.
  - Both papers use CMB dipole v = 369.82 ± 0.11 km/s toward (264.021°, 48.253°) (Planck 2018 I).
- Siewert et al. 2021 (A&A 653, A9) values were NOT verified here.

## 5. Axis of Evil
- Land & Magueijo 2005, PRL 95, 071301 (arXiv:astro-ph/0502237), WMAP1 maps: ℓ=2,3 planar multipoles aligned with **(b,l) ≈ (60°, −100°)**, i.e. (l,b) ≈ (260°, +60°); ℓ=5 preferred axis (b,l) = (50°, −91°) [= (269°, +50°)]. Alignment ℓ=2–5 rejects isotropy at >99.9%.
- Copi et al. and Planck 2018 VII axis values were NOT verified here (not fetched).

## Not included in data_small.tar.gz (>100 MB)
smica_2048 map (1.9 GB) and common mask (193 MB) stay on the box only.
