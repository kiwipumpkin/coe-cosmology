# Cycle 3c: noise-aware consistency test (fixed full-sample fits at each slice's own distances)

## Method
- **Models:** Chris's Newtonian integrator via the cycle3b tables, full-sample parameters held fixed with no refitting.
  - Ω_s = 0.049 ("omb"): D 49 Gpc.
  - Ω_s free ("free"): D 94 Gpc.
  - Each cluster gets the model A at its own longitude and distance (r = cz/71.6). These are averaged in the same sliding 90° cones as the data and normalised to the slice mean.
- **Data:** rebuilt Migkas 2020 sliding-cone A(ℓ) profile (5° steps, L–T slope fixed at the slice's own fit), normalised to its mean.
- **Errors:** 1000 cluster bootstraps. The overlapping cones are handled with the full bootstrap covariance.
  - The effective number of independent points, (tr C)²/tr(C²), is 4.6 / 5.2 / 5.3 / 5.2 (z<0.03 / 0.03–0.07 / z>0.07 / all), out of 72 cone positions.
- **Main χ²:** on the 1ℓ–3ℓ Fourier coefficients (6 numbers) with the bootstrap covariance.
  - Secondary χ²: on the profile using the leading eigenmodes of the bootstrap covariance (99 % of variance), Hartlap-corrected.
  - The isotropic flat line gets the same treatment.
- **Mocks:** 1000 per slice and model. For each mock, the L–T prediction is multiplied by the model A for each cluster and by a randomly permuted real residual. This keeps the real scatter and removes the real sky pattern. The profile is rebuilt with the identical pipeline.
  - Isotropic mocks (A = 1) are the null.
  - The mocks give calibrated χ² p-values and the noise distribution of the 2ℓ phase, the max/min longitudes and the 2ℓ amplitude.

## Results per slice
- **z<0.03 (n 32).**
  - Bootstrap-covariance χ²: omb 23.7 for 6, free 28.2 for 6, isotropic 20.8 for 6 (nominal p 6e-4 / 9e-5 / 2e-3).
  - With 32 clusters the Gaussian approximation fails. Mock-calibrated p: omb 0.058 (1.9σ), free 0.027 (2.2σ), isotropic 0.085 (1.7σ).
  - The noise spread of the 2ℓ phase is 45° rms. Observed 60° vs predicted 21°: a miss this large or larger happens in 41 % of mocks.
  - Observed 2ℓ amplitude 0.116 vs mock median 0.135: noise alone produces this (p 0.61).
  - **Verdict:** noise-dominated. Mild tension at about 2σ for every model, isotropy included. This slice can't discriminate.
- **0.03–0.07 (n 113).**
  - χ²: omb 13.1 for 6 (p 0.041; mock-calibrated 0.09, 1.7σ), free 11.8 (p 0.066; mock-calibrated 0.13, 1.5σ), isotropic 15.8 (p 0.015; mock-calibrated 0.046, 2.0σ).
  - The 154° phase, 41° from the predicted 15° (modulo 180), is matched or exceeded by noise in 24 % (omb) and 22 % (free) of mocks. Noise spread is about 36° rms.
  - The 2ℓ amplitude is the odd part: observed 0.151 vs model-mock median 0.070 (p 0.038, 2.1σ) and isotropic-mock median 0.061 (p 0.011).
  - **Verdict:** consistent with the model at 1.5–1.7σ. The phase shift is ordinary noise (about 1.2σ). The slice carries more 2ℓ power than the model predicts (about 2σ).
- **z>0.07 (n 168).**
  - χ²: omb 6.1 for 6 (p 0.41; mock-calibrated 0.45, 0.8σ), free 14.1 (p 0.029; mock-calibrated 0.05, 2.0σ), isotropic 4.9 (p 0.56; mock-calibrated 0.63).
  - Phase: observed 15° vs predicted 10°; misses this large happen in 80 % of mocks.
  - Amplitude: 0.064 vs mock median 0.055 (p 0.40).
  - **Verdict:** omb consistent (0.8σ). Free in 2σ tension. Isotropy is also fully consistent, so this slice cannot tell the model from no anisotropy.
- **All (n 313).**
  - χ²: omb 4.0 for 6 (p 0.68; mock-calibrated 0.69), free 6.7 (p 0.35; mock-calibrated 0.39), isotropic 12.7 (p 0.049; mock-calibrated 0.08, 1.8σ).
  - Phase: observed 8° vs predicted 12–14°; misses this large happen in 80–84 % of mocks.
  - Amplitude: 0.084 vs mock median 0.054 (p 0.13).
  - **Verdict:** both fits consistent (0.4σ and 0.9σ). The fixed omb model beats isotropy by Δχ² 8.7 with no parameters refitted to this profile, though these parameters were fitted to the Migkas full sample originally.
- **The three disjoint slices combined** (Fisher's method on the mock-calibrated p-values):
  - omb: p 0.059 (1.9σ).
  - free: p 0.008 (2.6σ).
  - isotropic: p 0.062 (1.9σ).
- **Max/min longitudes.**
  - The minimum is always consistent: mock p 0.25–0.99.
  - The maximum shows small p-values (0.000–0.07), but that is an artefact. The data and the model both have two peaks about 180° apart. The data's highest peak is the one near 180–210°, while the model's is near 25°. Noise swaps which peak is highest. The 2ℓ phase is the robust measure, and it is consistent everywhere.
- **Point-curve variant** (Chris's convention: unbinned model vs binned data), Fourier χ² for 6, omb / free:
  - z<0.03: 27.7 / 32.6.
  - 0.03–0.07: 14.8 / 13.4.
  - z>0.07: 6.8 / 17.7.
  - All: 3.1 / 7.3.
  - The conclusions are unchanged.

## Overall verdict
- **Correction to cycle 3b's wording.** Chris's points hold. Once the noise is modelled, R² was misleading: no slice is in more than about 2σ tension with the Ω_s = 0.049 fit, and the 154° phase in the middle slice is a roughly 1-in-4 noise fluctuation.
- **But the data don't discriminate much either.** Isotropy is equally consistent in z<0.03 and z>0.07, and only 1.8–2σ disfavoured overall.
- **Ω_s free fares slightly worse** (2σ at z>0.07; 2.6σ combined).
- **The one feature neither model reproduces** is the excess 2ℓ amplitude in the 0.03–0.07 slice (about 2σ).

## Files
- `noise_test.py` → `noise_test.json`: bootstraps, covariances, χ², mocks, phase/extremum/amplitude p-values.
- `make_fig.py` → `cycle3c_summary.png` (viewed):
  - top row: the model at each slice's distances inside the data's 16–84 % bootstrap band, with χ² and p in the legends;
  - bottom row: mock 2ℓ-phase distributions (model and isotropic) with the observed and predicted phases.

## Caveats
- Migkas 2020 L_X–T sample only. The profile rebuild uses a simplified fit (intrinsic scatter tuned per cone, E(z) from flat ΛCDM with Ω_m = 0.3).
- The mocks keep observed temperatures and permute log-residuals. That ignores the correlation between residuals and sky position caused by the real signal, and any latitude or absorption systematics.
- The bootstrap covariance is unreliable for n = 32. Use the mock-calibrated values there.
- The Fourier χ² keeps 1ℓ–3ℓ only. The eigenmode version (19–35 modes) gives the same ranking.
- The model is 2D (Galactic plane): clusters are assigned by longitude only.
- Parameters were fitted to the digitised Migkas full-sample profile. The "all" consistency is therefore not an independent test.
