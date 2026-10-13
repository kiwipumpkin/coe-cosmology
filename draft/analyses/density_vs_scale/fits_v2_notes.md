# v2 fits: |δ| = A·R^−β (R in Mpc), least squares in log-log (fits_v2.py → fits_v2.json)
The ± values are OLS standard errors. They are optimistic: the points mix tracers and definitions, and the 2M++ points are correlated.

| fit | n | β | scatter | 1 Gpc: peak ÷ fit | 2 Gpc | 4 Gpc |
|-|-|-|-|-|-|-|
| Over: all incl. R200 + 2M++ | 20 | 2.02 ± 0.15 | 0.29 dex | 36× | 119× | 520× |
| **Over: no R200 (chart line)** | 18 | **1.81 ± 0.22** | 0.28 | **19×** | **53×** | **200×** |
| Over: literature only, no R200 | 12 | 1.62 ± 0.38 | 0.33 | 8.5× | 21× | 70× |
| Over: literature + R200 | 14 | 2.00 ± 0.21 | 0.34 | 30× | 98× | 420× |
| Over: 2M++ alone | 6 | 1.86 ± 0.07 | 0.04 | 26× | 77× | 300× |
| Under: all incl. disputed + 2M++ | 17 | 0.53 ± 0.17 | 0.27 | 0.8× | 1.0× | 1.5× |
| **Under: no disputed (chart line)** | 12 | **0.66 ± 0.23** | 0.29 | **1.2×** | **1.6×** | **2.7×** |
| Under: literature only, no disputed | 6 | 0.54 ± 0.21 | 0.18 | 0.6× | 0.7× | 1.1× |
| Under: literature + disputed | 11 | 0.54 ± 0.15 | 0.18 | 0.65× | 0.76× | 1.2× |
| **Under: 2M++ alone (one survey, one method)** | 6 | **1.20 ± 0.24** | 0.16 | **9.8×** | **18×** | **45×** |
| Under, saturating form \|δ\|/(1−\|δ\|), no disputed | 12 | 1.06 ± 0.40 | 0.50 | 1.7× | 2.9× | 6.3× |

How fast ΛCDM's expected spread falls (σ ∝ R^−β, linear top-hat): β = 0.67 at 10 Mpc, 0.92 at 30, 1.25 at 100, 1.60 at 300, 1.78 at 1 Gpc, 1.83 at 2 Gpc, 1.86 at 4 Gpc. Pure white noise gives β = 1.5.
