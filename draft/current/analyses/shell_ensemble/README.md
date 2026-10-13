# Shell ensemble: random cluster+void web around the CoE

Scripts: `ensemble.py` (runs, `python3 ensemble.py 1000`), `counts.py` (counts/table), `figs.py` (figures). Raw per-realization arrays: `raw_*.npz`; summaries: `ensemble.json`, `counts.json`.

- Tidal code: `field()` from ../scan.py unchanged; `field3(alpha,beta,D)` is its 3D generalisation. beta=0 check: max relative difference 2.7e-15 over 200 random (alpha, D).
- kg conversion: KG = 7.7e51/0.0039 = 1.9744e+54 kg per unit k (same constant as make_figs.py, pos_free.py, web_floor.py, filament_demo.py). Equivalent to k = GM/(R_vis^2 c H0) with H0 = 70 (gives 7.74e51 kg at k = 0.0039).
- Geometry: CoE 58.2 Gpc from observer at ell = 195.16 deg (cmb-cancel convention: theta = ell - phi = 0 points away from O). Shell 43.9-72.5 Gpc around the CoE; positions uniform in shell volume, rejected if inside the 14.3 Gpc vision sphere. CoE at 15.2 deg instead gives the same statistics (best-scaled median 233.3 vs 233.8).
- Masses: |M| log-uniform within each level; half clusters, half voids, void deficits rescaled so net contrast = 0 exactly. Residual = R - sum k_i field3_i (scan.py sign convention). Best-scaled: one free multiplier s per realization (least squares, sign free), so it lowers RMS in 1000/1000 by construction.
- m3-5 and m1-2 'combined' = quadrature sum of A_m (A_m as in scan.py modes()).

Baseline (no web): RMS = 237.3× floor; |m3–5| = 322.6×, |m1–2| = 32.7× (quadrature sums of A_m). 1000 realizations per setting.

Mass scale → structure: radius of a δ=+1 cluster / δ=−0.5 void holding |M| (ρ̄ = 2.9e-27 kg/m³). Check: 7.6e51 kg → 2.77 / 3.49 Gpc.

|level|kg per object|cluster radius δ=+1 (Gpc)|void radius δ=−0.5 (Gpc)|
|-|-|-|-|
|(a) 1e48–1e49|1.0e+48–1.0e+49|0.14–0.30|0.18–0.38|
|(b) 3e49–3e50|3.2e+49–3.2e+50|0.45–0.96|0.56–1.21|
|(c) 1e51–1e52|1.0e+51–1.0e+52|1.41–3.04|1.78–3.83|

### objects with centre >14.3 Gpc from observer (spec)

|N|level|median RMS/floor|5th pct|min|best-scaled median|best 5th pct|best min|RMS lower / higher (placed)|m3–5 lower|m1–2 higher|3–5 lower & 1–2 not higher|best: 3–5 lower|best: 1–2 higher|best: both|
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|200|(a) 1e48–1e49|237.3|237.1|235.7|234.6|214.4|188.9|521 / 479|518|488|222|903|823|97|
|200|(b) 3e49–3e50|237.4|232.5|221.6|234.5|215.0|182.4|474 / 526|516|630|145|906|830|86|
|200|(c) 1e51–1e52|326.3|240.5|198.7|234.3|214.9|190.2|36 / 964|292|998|0|917|828|93|
|500|(a) 1e48–1e49|237.3|236.9|234.7|232.9|212.1|175.4|509 / 491|504|493|233|937|837|108|
|500|(b) 3e49–3e50|237.8|229.1|216.0|233.6|212.1|172.4|438 / 562|488|699|131|925|834|110|
|500|(c) 1e51–1e52|525.1|289.3|216.6|233.8|214.2|183.3|6 / 994|149|1000|0|936|826|122|

### whole object outside vision sphere (centre > 14.3 Gpc + own radius)

|N|level|median RMS/floor|5th pct|min|best-scaled median|best 5th pct|best min|RMS lower / higher (placed)|m3–5 lower|m1–2 higher|3–5 lower & 1–2 not higher|best: 3–5 lower|best: 1–2 higher|best: both|
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|200|(a) 1e48–1e49|237.3|237.1|235.9|234.2|214.2|174.8|509 / 491|496|485|198|920|828|97|
|200|(b) 3e49–3e50|237.3|232.8|220.3|234.0|214.9|140.2|454 / 546|500|599|145|926|855|75|
|200|(c) 1e51–1e52|289.3|234.4|194.0|233.7|212.7|173.4|77 / 923|311|1000|0|923|859|73|
|500|(a) 1e48–1e49|237.3|236.8|234.6|233.5|213.2|152.4|489 / 511|495|502|218|921|818|116|
|500|(b) 3e49–3e50|237.5|229.9|211.7|233.3|211.4|172.6|463 / 537|507|652|148|923|837|101|
|500|(c) 1e51–1e52|399.2|264.5|167.3|232.8|211.2|164.3|15 / 985|209|1000|0|934|855|90|


## E[RMS^2] = RMS0^2 + Var(web) check (floor^2 units)

|setting|mean RMS^2|RMS0^2 + mean Var(web)|median web RMS/floor|
|-|-|-|-|
|bare_N200_a|5.63e+04|5.63e+04|0.2|
|bare_N200_b|6.093e+04|6.08e+04|7.6|
|bare_N200_c|1.819e+06|1.82e+06|220.7|
|bare_N500_a|5.629e+04|5.63e+04|0.5|
|bare_N500_b|6.962e+04|6.963e+04|15.2|
|bare_N500_c|4.448e+06|4.449e+06|473.9|
|buf_N200_a|5.629e+04|5.629e+04|0.2|
|buf_N200_b|5.666e+04|5.669e+04|6.5|
|buf_N200_c|1.378e+05|1.383e+05|170.5|
|buf_N500_a|5.629e+04|5.629e+04|0.5|
|buf_N500_b|5.706e+04|5.725e+04|13.6|
|buf_N500_c|2.586e+05|2.575e+05|327.8|
|bare_N500_c_CoE_at_15.2|4.401e+07|4.401e+07|495.7|

Figures: `hist_rms_bare.png`, `hist_rms_buf.png`, `example_realization.png`.
