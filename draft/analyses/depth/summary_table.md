# Tidal influence depth ('multipole horizon')

Visible radius d = 14.3 Gpc. Leftover = RMS of the CMB residual after SR + central mass, in units of the 1e-5 floor (start 237.3×).

## A) Best leftover vs R_cut (objects allowed only within R_cut of us)

|series|16|18|20|22.5|25|28.6|32|36|42.9|50|57.2|final|flattens (≤5% above final)|90% of improvement by|95% by|
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|4 objects, masses or voids|179.7|130.6|104.8|87.0|78.4|73.7|73.0|73.0|73.0|73.0|73.0|73.0×|26.4 Gpc (1.84 d)|22.2 Gpc (1.55 d)|24.2 Gpc (1.69 d)|
|4 objects, masses only|179.7|130.6|104.8|87.0|78.4|73.7|73.0|73.0|73.0|||73.0×|26.4 Gpc (1.84 d)|22.2 Gpc (1.55 d)|24.2 Gpc (1.69 d)|
|3 objects, masses or voids|192.8|160.8|146.9|137.3|130.7|124.0|119.1|114.2|106.7|102.4|102.3|102.3×|42.3 Gpc (2.96 d)|34.7 Gpc (2.43 d)|40.8 Gpc (2.85 d)|
|4 objects, masses or voids (CoE interior allowed)|175.8|123.2|91.1|69.6|60.5|54.1|52.4|52.4|52.4|52.4|49.0|49.0×|52.0 Gpc (3.64 d)|23.0 Gpc (1.61 d)|26.2 Gpc (1.83 d)|
|4 objects, masses only (CoE interior allowed)|179.7|130.6|102.2|82.7|73.2|69.0|68.9|68.9|68.9|68.9|68.9|68.9×|25.7 Gpc (1.80 d)|22.1 Gpc (1.55 d)|23.9 Gpc (1.67 d)|
|3 objects, masses or voids (CoE interior allowed)|191.2|160.8|146.9|137.3|130.7|124.0|119.1|114.2|106.7|102.4|89.0|89.0×|54.8 Gpc (3.83 d)|47.6 Gpc (3.33 d)|53.2 Gpc (3.72 d)|

Series without '(CoE interior allowed)' keep every object ≥ D−d = 43.9 Gpc from the CoE (same region as the random web). Curves are made non-increasing (a larger R_cut can always reuse the smaller-R_cut solution).

## B1) Random web, k=0, N=200, medium, 50% voids (1,000 layouts, centres > 14.3 Gpc), split into distance shells

|shell (Gpc)|objects per layout|share of web effect (power)|share m=1–2|share m=3–5|pattern RMS per 1e51 kg (× floor)|m=1 / m=2 / m=3 / m=4 / m=5 amplitude per 1e51 kg (× floor)|median change in leftover from this shell alone (× floor)|
|-|-|-|-|-|-|-|-|
|14.3–20|14.8|99.69%|94.5%|97.9%|430.0|76.92 / 80.64 / 80.81 / 79.85 / 76.60|4.97|
|20–25|20.6|0.21%|3.4%|1.6%|16.9|9.37 / 15.28 / 11.15 / 7.83 / 5.79|1.95|
|25–30|30.0|0.06%|1.1%|0.3%|7.3|3.86 / 7.66 / 4.77 / 2.71 / 1.60|0.96|
|30–35|41.3|0.02%|0.5%|0.1%|4.1|1.96 / 4.71 / 2.32 / 1.19 / 0.57|0.56|
|35–42.9|93.4|0.02%|0.4%|0.1%|2.2|0.96 / 2.70 / 1.14 / 0.49 / 0.19|0.42|

## B2) Depth = radius containing 90% / 95% of the random web's effect (web extended to 400 Gpc, uniform in volume, CoE interior excluded)

|objects start at|all modes (leftover) R90 / R95|m=1–2 R90 / R95|m=3–5 R90 / R95|
|-|-|-|-|
|14.3 Gpc (centres anywhere outside; near-ring objects dominate)|14.8 / 15.3 Gpc (1.04 / 1.07 d)|18.5 / 21.8 Gpc (1.29 / 1.53 d)|16.8 / 18.3 Gpc (1.17 / 1.28 d)|
|15.0 Gpc (14.3 + typical object radius) ← headline|18.6 / 21.1 Gpc (1.30 / 1.47 d)|23.4 / 28.6 Gpc (1.64 / 2.00 d)|19.4 / 21.5 Gpc (1.36 / 1.50 d)|
|15.5 Gpc (larger voids)|20.8 / 24.1 Gpc (1.45 / 1.69 d)|25.6 / 31.5 Gpc (1.79 / 2.20 d)|20.6 / 22.9 Gpc (1.44 / 1.60 d)|
|B1 ensemble as run (14.3–42.9 Gpc)|14.7 / 15.1|17.7 / 20.3|16.4 / 17.8|

## C) Analytic falloff vs measured

|mode|expected amplitude falloff|measured slope 60–180 Gpc|measured slope 15–30 Gpc|far-field depth R90 / R95 (× inner edge r0)|
|-|-|-|-|-|
|m=1|(R_vis/r)^4|-4.01|-4.98|1.58 / 1.82|
|m=2|(R_vis/r)^3|-2.95|-3.95|2.15 / 2.71|
|m=3|(R_vis/r)^4|-3.93|-4.83|1.58 / 1.82|
|m=4|(R_vis/r)^5|-4.91|-5.76|1.39 / 1.53|
|m=5|(R_vis/r)^6|-5.88|-6.70|1.29 / 1.39|
